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
						   s::Int,
						   ω_::Int)

	# Remove item s-1
	ω_ += prob.W[1,seq[s-1]]
	sol.z -= prob.P[:,seq[s-1]]
	sol.X[seq[s-1]] = 0

	# Remove item s
	sol.z -= sol.X[seq[s]] * prob.P[:,seq[s]]
	sol.X[seq[s]] = 0

	if prob.W[1,seq[s]] <= ω_
		# Item s is inserted
		addItem!(prob, sol, seq[s])
		ω_ -= prob.W[1,seq[s]]

		if ω_ > 0
			# A fraction of item s-1 is inserted 
			addBreakItem!(prob, sol, ω_, seq[s-1])
		end

	else # Item s remains the break item
		addBreakItem!(prob, sol, ω_, seq[s])
		s = s-1 # The position of the break item changes
	end

	return sol, s, ω_
end

# The break item is swapped with an item that is not in the knapsack
function swapWithItemNotInBag(prob::_MOMKP,
						   	  seq::Vector{Int},
						   	  sol::Solution,
						   	  s::Int,
						   	  ω_::Int)

	# The item in position s is removed
	sol.z -= sol.X[seq[s]] * prob.P[:,seq[s]]
	sol.X[seq[s]] = 0

	if prob.W[1,seq[s+1]] <= ω_
		# Item s+1 can be inserted 
		addItem!(prob, sol, seq[s+1])
		ω_ -= prob.W[1,seq[s+1]]

		if ω_ > 0
			# A fraction of item s is inserted
			addBreakItem!(prob, sol, ω_, seq[s])
		end
		s = s+1 # The position of the break item changes

	else # Item s+1 becomes the break item
		addBreakItem!(prob, sol, ω_, seq[s+1])
	end

	return sol, s, ω_
end

# Computes the LP relaxation using the parametric method 
function parametricMethod(prob::_MOMKP,
						  L::Vector{Solution},
						  transpositions::Vector{Transposition},
						  seq::Vector{Int},
						  pos::Vector{Int})

	# Builds the initial solution
	sol, s, ω_ = buildSolution(prob, seq)

	upperBound = DualBoundSet{Float64}()
	updateBoundSets!(upperBound, L, 1//1, sol, seq[s])

	numberCasesIdenticalWeights = 0

	for iter in 1:length(transpositions)

		# Multiple identical critical weights
		if length(transpositions[iter].pairs) > 1

			numberCasesIdenticalWeights += 1

			# Positions corresponding to each transposition
			positions = [(min(pos[i], pos[j]), max(pos[i], pos[j])) 
						for (i,j) in transpositions[iter].pairs]
			sort!(positions)
			
			# Identification of the modified subsequences
			subsequences = identifySubsequences(positions)
			
			for (start,finish) in subsequences
					
				# The subsequence is reversed 
				seq[start:finish] = seq[finish:-1:start]
					
				if start <= s && finish >= s # The solution is modified
					sol, s, ω_ = reoptSolution(prob, seq, start, finish, sol, ω_)
					if ω_ > 0
						addBreakItem!(prob, sol, ω_, seq[s])
					end
				end
					
				updatePositions!(seq, pos, start, finish)
			end 

			updateBoundSets!(upperBound, L, transpositions[iter].λ, sol, seq[s]) 
		else

			(i,j) = transpositions[iter].pairs[1]
			k = min(pos[i], pos[j])

			if k == s-1   # Swap items s-1 and s

				sol, s, ω_ = swapWithItemInBag(prob, seq, sol, s, ω_)

			elseif k == s # Swap items s and s+1

				sol, s, ω_ = swapWithItemNotInBag(prob, seq, sol, s, ω_)
				
			end

			# Update the sequence and positions
			tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
			seq[pos[i]] = i ; seq[pos[j]] = j

			updateBoundSets!(upperBound, L, transpositions[iter].λ, sol, seq[s])
		end
	end

	#println("\tNumber of cases of identical critical weights : ", numberCasesIdenticalWeights)

	return upperBound
end
