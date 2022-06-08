################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Fonctions auxiliaires                                                #
################################################################################

include("dataStructures.jl")
include("orderedList.jl")

# ----- DOMINANCE ------------------------------------------------------------ #
# Returns true if x dominates y
# -- x and y are points 
function dominates(x, y, opt::Optimisation=MAX)
    if opt == MIN
        return ((x[1] <= y[1] && x[2] <  y[2])
             || (x[1] <  y[1] && x[2] <= y[2])
             || (x[1] == y[1] && x[2] == y[2])) # No duplicates
    else
        return ((x[1] >= y[1] && x[2] >  y[2])
             || (x[1] >  y[1] && x[2] >= y[2])
             || (x[1] == y[1] && x[2] == y[2])) # No duplicates
    end
end

# -- x and y are solutions 
function dominates(x::Solution, y::Solution, opt::Optimisation=MAX)
    if opt == MIN
        return ((x.z[1] <= y.z[1] && x.z[2] <  y.z[2])
             || (x.z[1] <  y.z[1] && x.z[2] <= y.z[2])
             || (x.z[1] == y.z[1] && x.z[2] == y.z[2])) # No duplicates
    else
        return ((x.z[1] >= y.z[1] && x.z[2] >  y.z[2])
             || (x.z[1] >  y.z[1] && x.z[2] >= y.z[2])
             || (x.z[1] == y.z[1] && x.z[2] == y.z[2])) # No duplicates
    end
end

# ----- INITIALISATION ------------------------------------------------------- #
# Computes the utilities (profit-to-weight ratios) for both objective functions
function utilities(prob::_MOMKP)

	n = size(prob.P)[2]

	r1 = [prob.P[1,j]//prob.W[1,j] for j in 1:n]
	r2 = [prob.P[2,j]//prob.W[1,j] for j in 1:n]
	return r1, r2
end

# Computes the critical weights
function criticalWeights(prob::_MOMKP,
						 r1::Vector{Rational{Int}},
						 r2::Vector{Rational{Int}})

	n       = size(prob.P)[2]
	weights = Rational{Int}[]
	pairs   = Tuple{Int,Int}[]

	nbTransp = 0

	# Computes the critical weight for each pair of items (i,j)
	for i in 1:n
		for j in i+1:n

			if !(r1[i] == r1[j] || r2[i] == r2[j]) 

				λ = (r2[j] - r2[i])//(r1[i] - r2[i] - r1[j] + r2[j])

				if λ > 0 && λ < 1
					nbTransp += 1 
					push!(weights, λ)
					push!(pairs, (i,j))
				end
			end
		end
	end

	#println("\tNumber of transpositions : ", nbTransp, " / ", Int(n*(n-1)/2))

	# Sorts the critical weights and associated item pairs in decreasing order
	perm = sortperm(weights, rev=true)
	return weights[perm], pairs[perm]
end

# Returns a list containing all the distinct λ values in decreasing order and
# the associated item pair(s)
function transpositionPreprocessing(weights::Vector{Rational{Int}},
									pairs::Vector{Tuple{Int,Int}})

	transpositions = Transposition[]
	nbIdenticalCriticalWeights = 0

	iter = 1
	while iter <= length(weights)

		# There are multiple identical critical weights λ
		if iter < length(weights) && weights[iter] == weights[iter+1]

			nbIdenticalCriticalWeights += 1 

			transp = Transposition(weights[iter])

			while iter < length(weights) && weights[iter] == weights[iter+1]
			push!(transp.pairs, pairs[iter])
			iter += 1
			end

			# The last occurence of λ
			push!(transp.pairs, pairs[iter])
			iter += 1

			push!(transpositions, transp)

		else
			push!(transpositions, Transposition(weights[iter], [pairs[iter]]))
			iter += 1
		end
	end

	#println("Number of identical critical weights : ", nbIdenticalCriticalWeights)
	return transpositions
end

# Computes the utilities, transpositions, initial sequence and positions 
function initialisation(prob::_MOMKP, method::Method)

    r1, r2 = utilities(prob)

	if method == DICHOTOMIC 

		return Initialisation(r1, r2, nothing, nothing, nothing)
	else
		# Item sequence 
		seq = sortperm(1000000*r1 + r2, rev=true) 

		if method == SIMPLEX 
			
			return Initialisation(r1, r2, nothing, seq, nothing)

		else # method == PARAMETRIC 

			# Item positions
			pos = sortperm(seq)          			  
			weights, pairs = criticalWeights(prob, r1, r2)
			transpositions = transpositionPreprocessing(weights, pairs)

			return Initialisation(r1, r2, transpositions, seq, pos)
		end 
	end 
end
# ----- SETTING VARIABLES ---------------------------------------------------- #
# Removes a variable from the sequence 
function removeFromSequence(seq::Vector{Int}, var::Int)

	newSeq = Vector{Int}(undef,length(seq)-1)

	# The variable is removed from the sequence
	inser = 1
	for i in 1:length(seq)
		if seq[i] != var
			newSeq[inser] = seq[i] 
			inser += 1 
		end
	end
	return newSeq 
end 

# Initial setvar 
function initialSetvar(prob::_MOMKP, 
					   init::Initialisation, 
					   method::Method)
	
	if method == PARAMETRIC_LP 
		transpInd = Vector{Int}[] 
		for i in 1:length(init.transpositions) 
			t = init.transpositions[i] 
			push!(transpInd, vcat([i], 1:length(t.pairs)))
		end 
		
		return SetVariables(Int[], Int[], transpInd, init.seq, init.pos)
	else
		return SetVariables(Int[], Int[], nothing, nothing, nothing)
	end 
end 

# Remove a variable from the initialisation data 
function setVariable(init::Initialisation,
					 parent_setvar::SetVariables, 
					 var::Int,					  # Variable to set 
					 val::Int, 					  # Value to set it to
					 method::Method)			  # Method for computing the UBS

	if method == PARAMETRIC_LP || method == PARAMETRIC_MT

		newTranspInd = Vector{Vector{Int}}() 
		newPos       = copy(parent_setvar.pos)

		# The variable is removed from the set of transpositions
		for i in 1:length(parent_setvar.transpInd)

			# Corresponding transposition 
			ind = parent_setvar.transpInd[i][1]
			t = init.transpositions[ind] 

			transpInd = Vector{Int}() 

			if length(t.pairs) > 1 

				swaps = Int[] 
				for indPair in parent_setvar.transpInd[i][2:end]
					if !(var in t.pairs[indPair])
						push!(swaps, indPair)
					end 
				end 

				if length(swaps) > 0 
					# The index of the critical weights and its associated 
					# pairs that do not contain var are inserted 
					push!(transpInd, ind)
					append!(transpInd, swaps)
					push!(newTranspInd, transpInd)
				end 
			else 

				if !(var in t.pairs[1])
					# The index of the critical weight and its only associated
					# pair are inserted 
					push!(transpInd, ind) 
					push!(transpInd, 1) 
					push!(newTranspInd, transpInd)
				end 
			end 
		end 

		# The variable is removed form the sequence 
		newSeq = removeFromSequence(parent_setvar.seq, var)

		# The positions of items after var in the sequence are diminished by 1
		for p in parent_setvar.pos[var]+1:length(parent_setvar.seq)
			newPos[parent_setvar.seq[p]] = 
				parent_setvar.pos[parent_setvar.seq[p]] - 1
		end

		if val == 0 
			return SetVariables(parent_setvar.setToOne, 
								vcat(parent_setvar.setToZero, [var]), 
								newTranspInd, newSeq, newPos)
		else 
			return SetVariables(vcat(parent_setvar.setToOne, [var]),
								parent_setvar.setToZero,
								newTranspInd, newSeq, newPos)
		end

	else # Dichotomic method and simplex algorithm
		if val == 0 
			return SetVariables(parent_setvar.setToOne, 
								vcat(parent_setvar.setToZero, [var]), 
								nothing, nothing, nothing)
		else 
			return SetVariables(vcat(parent_setvar.setToOne, [var]),
								parent_setvar.setToZero,
								nothing, nothing, nothing)
		end
	end 
end

# Remove a variable from the initialisation data (for the dichotomic method and 
# the simplex algorithm)
function setVariableInit(init::Initialisation, 
						 var::Int, 
						 method::Method)

	if method == SIMPLEX 

		n = length(init.r1) 

		newSeq = removeFromSequence(init.seq, var)

		r1 = Vector{Rational{Int}}(undef, n)
		r2 = Vector{Rational{Int}}(undef, n)

		for i in 1:n
			if i == var 
				r1[i] = -1
				r2[i] = -1 
			else 
				r1[i] = init.r1[i] 
				r2[i] = init.r2[i] 
			end 
		end 

		return Initialisation(r1, r2, nothing, newSeq, nothing)

	elseif method == DICHOTOMIC 

		n = length(init.r1) 

		r1 = Vector{Rational{Int}}(undef, n)
		r2 = Vector{Rational{Int}}(undef, n)

		for i in 1:n
			if i == var 
				r1[i] = -1
				r2[i] = -1 
			else 
				r1[i] = init.r1[i] 
				r2[i] = init.r2[i] 
			end 
		end 

		return Initialisation(r1, r2, nothing, nothing, nothing)
	end 
end 

# ----- SOLUTIONS ------------------------------------------------------------ #
# Add an item to a solution
function addItem!(prob::_MOMKP, 
				  sol::Solution, 
				  item::Int) 
	sol.X[item] = 1
	sol.z      += prob.P[:,item]
	sol.ω_     -= prob.W[1,item]
end

# Add a break item to a solution
function addBreakItem!(prob::_MOMKP,
					   sol::Solution,
					   item::Int)

	sol.X[item] = sol.ω_//prob.W[1,item]
	sol.z += sol.X[item] * prob.P[:,item]
end

# Computes the dantzig solution for a given sequence
function dantzigSolution(prob::_MOMKP, 
						 seq::Vector{Int}, 
						 setvar::Union{Nothing,SetVariables}=nothing)

	n   = size(prob.P)[2]
	sol = Solution{Float64}(prob)
	if !(setvar === nothing)
		for var in setvar.setToOne 
			addItem!(prob, sol, var)
		end 
	end 
	i   = 1

	while i <= length(seq) && prob.W[1,seq[i]] <= sol.ω_
		# L'objet est inséré
		addItem!(prob, sol, seq[i])
		i += 1
	end

	return sol, i
end

# Builds a solution including the break item
function buildSolution(prob::_MOMKP, 
					   seq::Vector{Int}, 
					   setvar::Union{Nothing,SetVariables}=nothing)

	n      = size(prob.P)[2]
	sol, s = dantzigSolution(prob, seq, setvar)

	if sol.ω_ > 0 && s <= length(seq)
		# Une fraction de l'objet s est insérée
		addBreakItem!(prob, sol, seq[s])
	end

	return sol, s
end

# Re-build part of a solution after a sequence reversal
function reoptSolution(prob::_MOMKP,
					   seq::Vector{Int},
					   start::Int,
					   finish::Int, 
					   sol::Solution) 

	# The items between start and finish in the sequence are removed 
	pos = finish 
	stop = false 
	while !stop && pos >= start

		item = seq[pos]
		
		if sol.X[item] < 1 && sol.X[item] > 0 
			# The break item is removed
			sol.z      -= sol.X[item] * prob.P[:,item]
			sol.X[item] = 0
			stop        = true
			
		elseif sol.X[item] == 1 
			# An item that was in the bag is removed 
			sol.z      -= prob.P[:,item]
			sol.ω_     += prob.W[1,item]
			sol.X[item] = 0
		end 
		pos -= 1 
	end

	# The items are inserted from start to finish 
	pos = start
	while prob.W[1,seq[pos]] <= sol.ω_

		# L'objet est inséré en entier
		addItem!(prob, sol, seq[pos])
		pos += 1
	end
	
	return sol, pos
end

# Returns true if the solution is integer 
function isInteger(sol::Solution, breakItem::Int)
	return (sol.X[breakItem] == 0//1 || sol.X[breakItem] == 1//1)
end

# Returns true if the solution is integer 
function isInteger(sol::Solution) 
	is_integer = true 
	for i in 1:length(sol.X)
		is_integer = is_integer && (sol.X[i] == 0//1 || sol.X[i] == 1//1)
	end 
	return is_integer
end

# Verifies the values for the objective function and the residual capacity
function verifySolution(prob::_MOMKP, sol::Solution) 
    objValue = [0,0]
    residualCapacity = prob.ω[1] 

    for i in 1:size(prob.P)[2] 
        if sol.X[i] == 1//1
            objValue += prob.P[:,i] 
            residualCapacity -= prob.W[1,i] 
        elseif sol.X[i] > 0 
            objValue += sol.X[i] * prob.P[:,i]
        end
    end

	@assert residualCapacity >= 0 "Solution non-admissible"
    println("z : ", sol.z == objValue)
    println("ω_ : ", sol.ω_ == residualCapacity)
end 

# ----- PROBLEMS ------------------------------------------------------------- #
# Groups together equivalent items 
function groupEquivalentItems(prob::_MOMKP)

    n      = size(prob.P)[2]
    r1, r2 = utilities(prob) 
    P1     = Int[] 
    P2     = Int[] 
    W      = Int[] 
    done   = Int[] # List of items that are equivalent to another item that has 
    # already been processed

    for i in 1:n       
        if !(i in done)

            p1 = prob.P[1,i]
            p2 = prob.P[2,i]
            w  = prob.W[1,i] 

            for j in i+1:n 

                if r1[i] == r1[j] && r2[i] == r2[j] 
                    # Items i and j are equivalent 
                    p1 += prob.P[1,j] 
                    p2 += prob.P[2,j] 
                    w  += prob.W[1,j] 
                    push!(done, j)
                end 
            end 

            push!(P1, p1)
            push!(P2, p2)
            push!(W, w) 
        end 
    end 

    P = Matrix(undef, 2, length(P1))
    for i in 1:length(P1) 
        P[1,i] = P1[i] ; P[2,i] = P2[i]  
    end 

    return _MOMKP(P, reshape(W, 1, length(W)), prob.ω)
end 

# TODO : tightness ratio 

# Transforms a multi-dimensional instance into a mono-dimensional instance by 
# only keep the first constraint 
function multiToMonoDimensional(prob::_MOMKP)
	m, n = size(prob.W)
	W = prob.W 
	for i in 1:m-1
		W = W[setdiff(1:end,2),:]
	end 
	return _MOMKP(prob.P, W, prob.ω[1:1])
end 

# ----- BOUND SETs ----------------------------------------------------------- #
# Returns true if sol is identical to the most recent solution in UB 
function identicalToPrevious(UB::Vector{Constraint}, 
							 sol::Solution) 
    return length(UB) > 0 && UB[end].point == sol.z 
end 

function identicalToPrevious(UB::DualBoundSet, y::Vector{Float64})
    return length(UB.points) > 0 && UB.points[end] == y
end 

# Updates the bound set by adding the point and corresponding constraint if it 
# is not already present and adding the solution to the list of integer 
# solutions if applicable 

# -- Parametric method for LP relaxation 
function updateBoundSets!(UB::DualBoundSet{Float64}, 
						  L::Vector{Solution{Float64}}, 
						  λ::Rational{Int}, 
						  sol::Solution{Float64}, 
						  breakItem::Int)
    
    if !identicalToPrevious(UB.constraints, sol)
        push!(UB.points, sol.z) 
		push!(UB.constraints, Constraint(λ, sol.z))
        if isInteger(sol, breakItem) 
            add!(L, Solution(sol.X[1:end], sol.z[1:end], sol.ω_)) 
        end 
    end 
end 

# -- Parametric method for Martello and Toth
function updateBoundSet!(UB::DualBoundSet{Float64}, 
						 λ::Union{Rational{Int},Float64}, 
						 y::Vector{Float64})

	if !identicalToPrevious(UB, y)
		add!(UB.points, y) 
		push!(UB.constraints, Constraint(λ, y))
	end 
end 

# -- Dichotomic method 
function updateBoundSets!(UB::DualBoundSet{Rational{Int}},  
						  L::Vector{Solution{Rational{Int}}},
						  λ::Rational{Int}, 
						  sol::Solution{Rational{Int}}, 
						  breakItem::Int)

	add!(UB.points, sol.z)
	push!(UB.constraints, Constraint(λ, sol.z))
	if isInteger(sol, breakItem) 
		add!(L, Solution(copy(sol.X), copy(sol.z), sol.ω_)) 
	end 
end

# -- Simplex algorithm 
function updateBoundSets!(UB::DualBoundSet{Float64},  
						  L::Vector{Solution{Float64}},
						  sol::Solution{Float64}, 
						  breakItem::Int)

	if !identicalToPrevious(UB, sol.z)
		push!(UB.points, sol.z) 

		# Constraint 
		if length(UB.constraints) == 0 
			@assert length(UB.points) == 1 "There are no constraints because the 
			first point has just been added"
			push!(UB.constraints, Constraint(1//1, sol.z))
		else 
			λ1 = abs(sol.z[2] - UB.points[end][2])
			λ2 = abs(sol.z[1] - UB.points[end][1])
			λ  = λ1/(λ1 + λ2)
			push!(UB.constraints, Constraint(λ, sol.z))
		end 

		if isInteger(sol, breakItem) 
			add!(L, Solution(copy(sol.X), copy(sol.z), sol.ω_)) 
		end 
	end 
end
