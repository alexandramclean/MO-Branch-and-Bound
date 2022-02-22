################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Fonctions auxiliaires                                                #
################################################################################

include("dataStructures.jl")

# ----- DOMINANCE ------------------------------------------------------------ #
# Returns true if x dominates y
function domine(x, y, opt::Optimisation=Max)
    if opt == Min
        return ((x[1] <= y[1] && x[2] < y[2])
            || (x[1] < y[1] && x[2] <= y[2])
            || (x[1] == y[1] && x[2] == y[2])) # No duplicates
    else
        return ((x[1] >= y[1] && x[2] > y[2])
            || (x[1] > y[1] && x[2] >= y[2])
            || (x[1] == y[1] && x[2] == y[2])) # No duplicates
    end
end

# ----- RATIOS AND CRITICAL WEIGHTS ------------------------------------------ #
# Computes the ratios for both objective functions
function ratios(prob::_MOMKP)

	n = size(prob.P)[2]
	r = Vector{Tuple{Rational{Int}, Rational{Int}}}(undef, n)

	r1 = [prob.P[1,i]//prob.W[1,i] for i in 1:n]
	r2 = [prob.P[2,i]//prob.W[1,i] for i in 1:n]

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

				λ = (r2[j] - r2[i])//(r1[i]-r2[i]-r1[j]+r2[j])

				if λ > 0 && λ < 1
					nbTransp += 1 
					push!(weights, λ)
					push!(pairs, (i,j))
				end
			end
		end
	end

	println("\tNumber of transpositions : ", nbTransp, " / ", n*(n-1)/2)

	# Sorts the critical weights and associated item pairs in decreasing order
	perm    = sortperm(weights, rev=true)
	return weights[perm], pairs[perm]
end

# Returns a list containing all the distinct λ values in decreasing order and
# the associated item pair(s)
function transpositionPreprocessing(weights::Vector{Rational{Int}},
					   				pairs::Vector{Tuple{Int,Int}})

	transpositions = Transposition[]

	iter = 1
	while iter <= length(weights)

		# There are multiple identical critical weights λ
		if iter < length(weights) && weights[iter] == weights[iter+1]

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
	return transpositions
end

# Updates the positions after reversing the subsequence between start and finish
function updatePositions!(seq::Vector{Int}, 
						  pos::Vector{Int}, 
						  start::Int, 
						  finish::Int) 

	if finish - start >= 1
		# Swaps the positions of the items at the edge of the subsequence 
		tmp              = pos[seq[start]] 
		pos[seq[start]]  = pos[seq[finish]]
		pos[seq[finish]] = tmp 
		
		# Appel récursif sur la sous-séquence privée des extrémités
		updatePositions!(seq, pos, start+1, finish-1)
	end
end

# ----- INITIALISATION ------------------------------------------------------- #
function initialisation(prob::_MOMKP)

	# Ratios and critical weights
	r1, r2 = ratios(prob)
	weights, pairs = criticalWeights(prob, r1, r2)
	
	# Identical critical weights are grouped together 
	transpositions = transpositionPreprocessing(weights, pairs)
	
	# Sorts the ratios in lexicographically decreasing order according to (r1, r2) 
	seq = sortperm(1000000*r1 + r2, rev=true) # Item sequence 
	pos = sortperm(seq)          			  # Item positions
	
	return transpositions, seq, pos
end

# ----- SOLUTIONS ------------------------------------------------------------ #
# Add an item to a solution
function addItem!(prob::_MOMKP, sol::Solution, item::Int)
	sol.X[item] = 1
	sol.z += prob.P[:,item]
end

# Add a break item to a solution
function addBreakItem!(prob::_MOMKP,
					   sol::Solution,
					   ω_::Int,
					   item::Int)

	sol.X[item] = ω_//prob.W[1,item]
	sol.z += sol.X[item] * prob.P[:,item]
end

# Computes the dantzig solution for a given sequence
function dantzigSolution(prob::_MOMKP, sequence::Vector{Int})

	n   = size(prob.P)[2]
	ω_  = prob.ω[1]
	sol = Solution(n)
	i   = 1

	while i <= n && prob.W[1,sequence[i]] <= ω_
		item = sequence[i]
		# L'objet est inséré
		addItem!(prob, sol, item)
		ω_ -= prob.W[1,item]
		i += 1
	end

	return sol, i, ω_
end

# Builds a solution including the break item
function buildSolution(prob::_MOMKP, seq::Vector{Int})

	n          = size(prob.P)[2]
	sol, s, ω_ = dantzigSolution(prob, seq)

	if ω_ > 0
		# Une fraction de l'objet s est insérée
		addBreakItem!(prob, sol, ω_, seq[s])
	end

	return sol, s, ω_
end

# Re-build part of a solution after a sequence reversal
function reoptSolution(prob::_MOMKP,
					   seq::Vector{Int},
					   start::Int,
					   finish::Int, 
					   sol::Solution,
					   s::Int,
					   ω_::Int)

	# Les objets entre start et fin dans la séquence sont retirés
	for pos in start:finish

		item = seq[pos]
		
		if sol.X[item] < 1 && sol.X[item] > 0 
			# L'objet cassé est retiré
			sol.z      -= sol.X[item] * prob.P[:,item]
			sol.X[item] = 0
			
		elseif sol.X[item] == 1 
			# Un objet inséré dans le sac est retiré
			sol.z      -= prob.P[:,item]
			ω_         += prob.W[1,item]
			sol.X[item] = 0
		end 
	end

	# Les objets sont insérés en partant de start dans la nouvelle séquence
	pos = start
	while prob.W[1,seq[pos]] <= ω_

		# L'objet est inséré en entier
		addItem!(prob, sol, seq[pos])
		ω_ -= prob.W[1,seq[pos]]

		pos += 1
	end

	return sol, pos, ω_
end
