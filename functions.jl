################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Fonctions auxiliaires                                                #
################################################################################

include("dataStructures.jl")

# ----- DOMINANCE ------------------------------------------------------------ #
# Retourne vrai si x domine y
function domine(x, y, opt::Optimisation=Max)
    if opt == Min
        return ((x[1] <= y[1] && x[2] < y[2])
            || (x[1] < y[1] && x[2] <= y[2])
            || (x[1] == y[1] && x[2] == y[2])) # Pas de doublons
    else
        return ((x[1] >= y[1] && x[2] > y[2])
            || (x[1] > y[1] && x[2] >= y[2])
            || (x[1] == y[1] && x[2] == y[2])) # Pas de doublons
    end
end

# ----- RATIOS ET POIDS CRITIQUES -------------------------------------------- #
# Calcul des ratios pour les deux fonctions objectif
function ratios(prob::_MOMKP)

	n  = size(prob.P)[2]
	r1 = Vector{Rational{Int}}(undef,n) 
	r2 = Vector{Rational{Int}}(undef,n)
	
	for i in 1:n 
		@assert prob.W[1,i] != 0 "Un objet ne peux pas avoir un poids de 0"
		r1 = [prob.P[1,i]//prob.W[1,i] for i in 1:n]
		r2 = [prob.P[2,i]//prob.W[1,i] for i in 1:n]
	end
	
	return r1, r2
end

function lambda(r1, r2, i, j) 
	return r2[j] - r2[i], r1[i] - r2[i] - r1[j] + r2[j] 
end 

# Calcul des poids critiques
function criticalWeights(prob::_MOMKP,
						 r1::Vector{Rational{Int}},
						 r2::Vector{Rational{Int}})

	n       = size(prob.P)[2]
	weights = Rational{Int}[]
	pairs   = Tuple{Int,Int}[]
	
	# Calcul des poids critiques pour chaque couple d'objets (i,j)
	for i in 1:n
		for j in i+1:n

			if (r1[i]-r2[i]-r1[j]+r2[j]) != 0 && !(r1[i] == r1[j] || r2[i] == r2[j])
				
				λ = (r2[j] - r2[i])//(r1[i]-r2[i]-r1[j]+r2[j])
				
				if λ > 0 && λ < 1
					push!(weights, λ)
					push!(pairs, (i,j))
				end
			end
		end
	end

	# Tri des poids critiques dans l'ordre décroissant
	perm = sortperm(weights, rev=true)
	return weights[perm], pairs[perm]
end

# Retourne une liste contenant une suite de valeurs de λ distinctes ainsi que 
# les transpositions correspondantes
function transpositionPreprocessing(weights::Vector{Rational{Int}}, 
					   				pairs::Vector{Tuple{Int,Int}})
	transpositions = Transposition[]
	iter = 1
	while iter <= length(weights) 
		# Il y a plusieurs poids critiques identiques 
		if iter < length(weights) && weights[iter] == weights[iter+1] 
			transp = Transposition(weights[iter]) 
			
			while iter < length(weights) && weights[iter] == weights[iter+1] 
				push!(transp.pairs, pairs[iter])
				iter += 1
			end
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
# Ajout d'un objet entier à une solution
function addItem!(prob::_MOMKP, sol::Solution, item::Int)
	sol.X[item] = 1
	sol.z += prob.P[:,item]
end

# Ajout d'un objet cassé à une solution
function addBreakItem!(prob::_MOMKP,
					   sol::Solution,
					   residualCapacity::Union{Int,Rational{Int}},
					   item::Int)

	sol.X[item] = residualCapacity//prob.W[1,item]
	sol.z += sol.X[item] * prob.P[:,item]
end

# Calcul de la solution dantzig pour une séquence donnée
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

# Construction d'une solution avec objet cassé
function buildSolution(prob::_MOMKP, sequence::Vector{Int})

	n = size(prob.P)[2]
	sol, s, residualCapacity = dantzigSolution(prob, sequence)

	if residualCapacity > 0
		# Une fraction de l'objet s est insérée
		addBreakItem!(prob, sol, residualCapacity, sequence[s])
	end

	return sol, s, residualCapacity
end

#= Exemple 
weights = [9//10, 85//100, 8//10, 8//10, 66//100, 65//100, 4//10, 4//10, 4//10]
pairs = [(1,5), (2,4), (2,5), (3,6), (4,5), (1,2), (1,4), (5,6), (3,5)]
@assert length(weights) == length(pairs) "Il doit y avoir autant de poids que de paires"
transpositions = transpositionPreprocessing(weights, pairs)
for t in transpositions
	println(t) 
end =#
