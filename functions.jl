################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Fonctions auxiliaires                                                #
################################################################################

# Calcul des ratios pour les deux fonctions objectif
function ratios(prob::_MOMKP)

	n  = size(prob.P)[2]
	r1 = Vector{Float64}(undef, n)
	r2 = Vector{Float64}(undef, n)
	for i in 1:n 
		r1[i] = prob.P[1,i]/prob.W[1,i] 
		r2[i] = prob.P[2,i]/prob.W[1,i]
	end
	return r1, r2
end

# Calcul des poids critiques
function criticalWeights(prob::_MOMKP, r1, r2)

	n       = size(prob.P)[2]
	weights = Float64[]
	pairs   = Tuple{Int64,Int64}[] 
	# Calcul des poids critiques pour chaque couple d'objets (i,j)
	for i in 1:n
		for j in i:n
			if r1[i]-r2[i]-r1[j]+r2[j] != 0 && !(r1[i] == r1[j] || r2[i] == r2[j])
				λ = (r2[j] - r2[i])/(r1[i]-r2[i]-r1[j]+r2[j])
				if λ > 0 && λ < 1
					push!(weights, λ)
					push!(pairs, (i,j))
				end
			end
		end
	end
	
	@assert length(weights) == length(pairs) "Il doit y avoir autant de paires 
	que de valeurs de λ"
	
	# Tri des poids critiques dans l'ordre décroissant
	perm = sortperm(weights, rev=true)
	return weights[perm], pairs[perm]
end

# Ajout d'un objet entier à une solution
function addItem!(prob::_MOMKP, sol::Solution, item) 
	sol.X[item] = 1
	sol.z += prob.P[:,item] 
end

# Ajout d'un objet cassé à une solution 
function addBreakItem!(prob::_MOMKP, sol::Solution, residualCapacity, item) 
	sol.X[item] = residualCapacity/prob.W[1,item] 
	sol.z += sol.X[item] * prob.P[:,item] 
end

# Calcul de la solution dantzig pour une séquence donnée
function dantzigSolution(prob::_MOMKP, sequence)

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

# Génère une fonction objectif pondérée
function weightedSum(prob::_MOMKP, λ)
    
    n = size(P.P)[2]
    # Construire la fonction objectif pondérée
    weightedObj = Vector{Float64}(undef,n)
    for j in 1:n
        weightedObj[j] = λ*prob.P[1,j] + (1 - λ)*prob.P[2,j]
    end
    return weightedObj
end

