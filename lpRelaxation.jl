################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Calcul paramétrique de la relaxation continue                        #
################################################################################

include("parametricMethodFunctions.jl")

# The break item is swapped with an item in the knapsack
function swapWithItemInBag(prob::_MOMKP,
						   seq::Vector{Int},
						   sol::Solution,
						   s::Int)

	# Remove item s-1
	sol.ω_ += prob.W[1,seq[s-1]]
	sol.z  -= prob.P[:,seq[s-1]]
	sol.X[seq[s-1]] = 0

	# Remove item s
	sol.z -= sol.X[seq[s]] * prob.P[:,seq[s]]
	sol.X[seq[s]] = 0

	if prob.W[1,seq[s]] <= sol.ω_
		# Item s is inserted
		addItem!(prob, sol, seq[s])

		if sol.ω_ > 0
			# A fraction of item s-1 is inserted 
			addBreakItem!(prob, sol, seq[s-1])
		end

	else # Item s remains the break item
		addBreakItem!(prob, sol, seq[s])
		s = s-1 # The position of the break item changes
	end

	return sol, s
end

# The break item is swapped with an item that is not in the knapsack
function swapWithItemNotInBag(prob::_MOMKP,
						   	  seq::Vector{Int},
						   	  sol::Solution,
						   	  s::Int)

	# The item in position s is removed
	sol.z -= sol.X[seq[s]] * prob.P[:,seq[s]]
	sol.X[seq[s]] = 0

	if prob.W[1,seq[s+1]] <= sol.ω_
		# Item s+1 can be inserted 
		addItem!(prob, sol, seq[s+1])

		if sol.ω_ > 0
			# A fraction of item s is inserted
			addBreakItem!(prob, sol, seq[s])
		end
		s = s+1 # The position of the break item changes

	else # Item s+1 becomes the break item
		addBreakItem!(prob, sol, seq[s+1])
	end

	return sol, s
end

# There are multiple identical critical weights 
function identicalCriticalWeights(prob::_MOMKP,
								  seq::Vector{Int},
								  pos::Vector{Int},
								  pairs::Vector{Tuple{Int,Int}},
								  sol::Solution{Float64},
								  s::Int)

	# Positions corresponding to each transposition
	positions = [(min(pos[i], pos[j]), max(pos[i], pos[j])) 
		for (i,j) in pairs]
	sort!(positions)

	# Identification of the modified subsequences
	subsequences = identifySubsequences(positions)

	for (start,finish) in subsequences

		# The subsequence is reversed 
		seq[start:finish] = seq[finish:-1:start]

		if start <= s && finish >= s # The solution is modified
			sol, s = reoptSolution(prob, seq, start, finish, sol)
			if sol.ω_ > 0
				addBreakItem!(prob, sol, seq[s])
			end
		end

		updatePositions!(seq, pos, start, finish)
	end 
	return sol, s 
end 

# ----- WITHOUT INTERRUPTION ------------------------------------------------- #
# Computes the LP relaxation using the parametric method 
function parametricMethod(prob::_MOMKP,           # Bi01KP instance
						  init::Initialisation,   # seq, pos, transpositions
						  setvar::SetVariables,   # 
						 ) where T<:Real

	# Creates copies of the sequence and positions as they will be modified 
	seq = setvar.seq[1:end] 
	pos = setvar.pos[1:end] 
	
	# Builds the initial solution
	sol, s = buildSolution(prob, seq, setvar)

	# The upper bound set is stored 
	UB = DualBoundSet{Float64}()
	# Stores the integer solutions found during the computation of the UBS 
	Lη = Vector{Solution{Float64}}()

	if s <= length(seq)
		updateBoundSets!(UB, Lη, 1//1, sol, seq[s])
	elseif length(seq) > 0 
		updateBoundSets!(UB, Lη, 1//1, sol, seq[s-1])
	end 

	numberCasesIdenticalWeights = 0

	iter = 1 
	while iter <= length(setvar.transpInd)

		# The first integer in the list is the index of the critical weight 
		# The others are the indices of the corresponding pairs 
		ind = setvar.transpInd[iter][1]

		λ     = init.transpositions[ind].λ 
		pairs = init.transpositions[ind].pairs[setvar.transpInd[iter][2:end]]

		# Multiple identical critical weights
		if length(pairs) > 1

			numberCasesIdenticalWeights += 1

			sol, s = identicalCriticalWeights(prob, seq, pos, pairs, sol, s)			
		else

			(i,j) = pairs[1]
			k = min(pos[i], pos[j])

			if k == s-1   # Swap items s-1 and s

				sol, s = swapWithItemInBag(prob, seq, sol, s)

			elseif k == s # Swap items s and s+1

				sol, s = swapWithItemNotInBag(prob, seq, sol, s)
				
			end

			# Update the sequence and positions
			tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
			seq[pos[i]] = i ; seq[pos[j]] = j
		end

		if s <= length(seq)
			updateBoundSets!(UB, Lη, λ, sol, seq[s]) 
		else 
			updateBoundSets!(UB, Lη, λ, sol, seq[s-1]) 
		end 

		iter += 1 
	end

	# Add the last constraint : z2 <= sol.z[2]
	push!(UB.constraints, Constraint(0//1, sol.z))

	#println("\tNumber of cases of identical critical weights : ", numberCasesIdenticalWeights)
	return UB, Lη
end

# ----- WITH INTERRUPTION ---------------------------------------------------- #
# Tests all the shifted local nadir points after a new constraint has been added
# to determine whether the computation of the upper bound set can be interrupted 
# Updates toBeTested by removing the local nadir points that do not verify 
# the new constraint  
function testLocalNadirPoints(L::Vector{Solution{Float64}},
							  # The critical weight and points that define 
							  # the new constraint 
							  λ::Rational{Int},
							  a1::Vector{Float64},
							  a2::Vector{Float64},
							  # Indices of nadir points to be tested
							  toBeTested::Vector{Int},
							  numberNadirsLeft::Int, 
							  # The computation of the upper bound set can be 
							  # interrupted 
							  interrupt::Bool)

	is_interrupted = false 

	for i in 1:length(toBeTested)
		j = toBeTested[i] 
		if j != 0
			# Does the local nadir point verify the new constraint ?
			nadir = [L[j-1].z[1]+1., L[j].z[2]+1.]

			if λ*nadir[1] + (1-λ)*nadir[2] <= λ*a1[1] + (1-λ)*a1[2] 
				# Can the computation be interrupted ? 
				if interrupt && nadir[1] <= a2[1] && nadir[2] <= a1[2] 
					is_interrupted = true 
				end 
			else 
				toBeTested[i] = 0 
				numberNadirsLeft -= 1 
			end 
		end 
	end 
	return is_interrupted, numberNadirsLeft
end 

# Computes the LP relaxation using the parametric method 
# The upper bound set is not stored, the status of the corresponding node is 
# returned indicating whether it is dominated by the lower bound set L, optimal,
# or not pruned  
# If the computation of the upper bound set is interrupted then it is not 
# dominated by the lower bound set 
function parametricLPrelaxation(prob::_MOMKP,           # Bi01KP instance
						  	    L::Vector{Solution{T}}, # Lower bound set 
						  	    init::Initialisation,   # seq, pos, transpositions
						  	    setvar::SetVariables,   # 
						  	    interrupt::Bool = true  # The computation of the 
						  	    # upper bound set can be interrupted
						 	   ) where T<:Real

	# Creates copies of the sequence and positions as they will be modified 
	seq::Vector{Int} = setvar.seq[1:end] 
	pos::Vector{Int} = setvar.pos[1:end] 

	# Stores the integer solutions found during the computation of the UBS 
	Lη = Vector{Solution{Float64}}() 

	# Builds the initial solution
	sol, s = buildSolution(prob, seq, setvar)

	# If the solution is integer it is stored 
	if s <= length(seq)
		if isInteger(sol, seq[s])
			add!(Lη, Solution(sol.X[1:end], sol.z[1:end], sol.ω_))
		end 
	elseif length(seq) > 0 
		if isInteger(sol, seq[s-1])
			add!(Lη, Solution(sol.X[1:end], sol.z[1:end], sol.ω_))
		end 
	end

	status = NOTPRUNED

	# First point 
	firstPoint = copy(sol.z)

	# Used to determine if the computation of the upper bound set is interrupted
	a2             = sol.z
	toBeTested     = [i for i in 2:length(L)] 
	nbNadirsLeft   = length(toBeTested)
	is_interrupted = false 

	iter = 1 
	while !is_interrupted && iter <= length(setvar.transpInd)

		# The first integer in the list is the index of the critical weight 
		# The others are the indices of the corresponding pairs 
		ind = setvar.transpInd[iter][1]

		λ     = init.transpositions[ind].λ 
		pairs = init.transpositions[ind].pairs[setvar.transpInd[iter][2:end]]

		# Multiple identical critical weights
		if length(pairs) > 1

			sol, s = identicalCriticalWeights(prob, seq, pos, pairs, sol, s)			
		else

			(i,j) = pairs[1]
			k = min(pos[i], pos[j])

			if k == s-1   # Swap items s-1 and s

				sol, s = swapWithItemInBag(prob, seq, sol, s)

			elseif k == s # Swap items s and s+1

				sol, s = swapWithItemNotInBag(prob, seq, sol, s)

			end

			# Update the sequence and positions
			tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
			seq[pos[i]] = i ; seq[pos[j]] = j
		end

		# If the solution is integer it is stored 
		if s <= length(seq)
			if isInteger(sol, seq[s])
				add!(Lη, Solution(sol.X[1:end], sol.z[1:end], sol.ω_))
			end 

		elseif length(seq) > 0 
			if isInteger(sol, seq[s-1])
				add!(Lη, Solution(sol.X[1:end], sol.z[1:end], sol.ω_))
			end 
		end

		# The new constraint is defined by λ, a1 and a2 
		if length(L) > 0 
			a1 = sol.z
			if a1 != a2 
				is_interrupted, nbNadirsLeft = testLocalNadirPoints(L, λ, 
					a1, a2, toBeTested, nbNadirsLeft, interrupt)
				a2 = a1 
			end 
		end 

		iter += 1 
	end

	if !interrupt && length(L) > 0
		# Test the last constraint : z2 <= sol.z[2]
		is_interrupted, nbNadirsLeft = testLocalNadirPoints(L, 0//1, sol.z, a2, 
			toBeTested, nbNadirsLeft, interrupt)
	end 

	# Last point 
	lastPoint = sol.z

	# The upper bound set is dominated if there are no shifted local nadir 
	# points remaining that verify all the constraints 
	if length(L) > 0
		if nbNadirsLeft == 0
			status = DOMINANCE
		elseif length(Lη) == 1 && firstPoint == lastPoint == Lη[1].z 
			status = OPTIMALITY
		end 
	end 

	return Lη, status 
end

