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

	println("\tNombre de transpositions : ", nbTransp, " / ", n*(n-1)/2)

	# Sorts the critical weights and associated item pairs in decreasing order
	perm = sortperm(weights, rev=true)
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

# ----- SOLUTIONS ------------------------------------------------------------ #
# Add an item to a solution
function addItem!(prob::_MOMKP, sol::Solution, item::Int)
	sol.X[item] = 1
	sol.z += prob.P[:,item]
end

# Add a break item to a solution
function addBreakItem!(prob::_MOMKP,
					   sol::Solution,
					   residualCapacity::Int,
					   item::Int)

	sol.X[item] = residualCapacity//prob.W[1,item]
	sol.z += sol.X[item] * prob.P[:,item]
end

# Computes the dantzig solution for a given sequence
function dantzigSolution(prob::_MOMKP, sequence::Vector{Int})

	n                = size(prob.P)[2]
	residualCapacity = prob.ω[1]
	sol              = Solution(n)
	i                = 1

	while i <= n && prob.W[1,sequence[i]] <= residualCapacity
		item = sequence[i]
		# L'objet est inséré
		addItem!(prob, sol, item)
		residualCapacity -= prob.W[1,item]
		i += 1
	end

	return sol, i, residualCapacity
end

# Builds a solution including the break item
function buildSolution(prob::_MOMKP, sequence::Vector{Int})

	n = size(prob.P)[2]
	sol, s, residualCapacity = dantzigSolution(prob, sequence)

	if residualCapacity > 0
		# Une fraction de l'objet s est insérée
		addBreakItem!(prob, sol, residualCapacity, sequence[s])
	end

	return sol, s, residualCapacity
end

# Re-build part of a solution
function reoptSolution(prob::_MOMKP,
					   prev::Vector{Int},
					   seq::Vector{Int},
					   deb::Int,
					   sol::Solution,
					   s::Int,
					   residualCapacity::Int)

	# Les objets entre deb et s dans la séquence précédente sont retirés
	for pos in deb:s-1

		item = prev[pos]
		sol.X[item] = 0
		sol.z -= prob.P[:,item]
		residualCapacity += prob.W[1,item]

	end

	# L'objet cassé est retiré
	sol.z -= sol.X[prev[s]] * prob.P[:,prev[s]]
	sol.X[prev[s]] = 0

	# Les objets sont insérés en partant de deb dans la nouvelle séquence
	pos = deb
	while prob.W[1,seq[pos]] <= residualCapacity

		# L'objet est inséré en entier
		addItem!(prob, sol, seq[pos])
		residualCapacity -= prob.W[1,seq[pos]]

		pos += 1
	end

	if residualCapacity > 0
		addBreakItem!(prob, sol, residualCapacity, seq[pos])
	end

	return sol, pos, residualCapacity

end
