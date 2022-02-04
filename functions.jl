################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Fonctions auxiliaires                                                #
################################################################################

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
	r1 = Vector{Rational{Int}}(undef, n)
	r2 = Vector{Rational{Int}}(undef, n)
	for i in 1:n 
		r1[i] = prob.P[1,i]//prob.W[1,i] 
		r2[i] = prob.P[2,i]//prob.W[1,i]
	end
	return r1, r2
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
		for j in i:n
			if r1[i]-r2[i]-r1[j]+r2[j] != 0 && !(r1[i] == r1[j] || r2[i] == r2[j])
				λ = (r2[j] - r2[i])//(r1[i]-r2[i]-r1[j]+r2[j])
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

# ----- SOLUTIONS ------------------------------------------------------------ #
# Ajout d'un objet entier à une solution
function addItem!(prob::_MOMKP, sol::Solution, item::Int) 
	sol.X[item] = 1
	sol.z += prob.P[:,item] 
end

# Ajout d'un objet cassé à une solution 
function addBreakItem!(prob::_MOMKP, 
					   sol::Solution, 
					   residualCapacity::Rational{Int}, 
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

	n                = size(prob.P)[2]
	residualCapacity = prob.ω[1]
	sol              = Solution(n)
	i                = 1
	
	while residualCapacity > 0 && i <= n
		item = sequence[i]
		
		if prob.W[1,item] <= residualCapacity # L'objet est inséré en entier
			addItem!(prob, sol, item)
			residualCapacity -= prob.W[1,item]
			
		else # Une fraction de l'objet est insérée
			addBreakItem!(prob, sol, residualCapacity, item)
		end
		i += 1
	end
	
	# Position de l'objet cassé 
	if sol.X[sequence[i-1]] < 1
		s = i-1 
	else 
		s = i
	end
	return sol, s, residualCapacity
end

