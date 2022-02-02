################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Calcul paramétrique de la borne de Martello et Toth                  #
################################################################################

include("dataStructures.jl")
include("functions.jl")
include("listeOrdonnee.jl") 

# Borne de Martello et Toth
function uMT(prob::_MOMKP, 
			 seq,
	     	 sol, # Solution de Dantzig
	     	 residualCapacity, 
	     	 s) # Position de l'objet cassé 
	
	u0 = sol.z + [0,0] 
	u0[1] += residualCapacity * (prob.P[1,seq[s+1]]/prob.W[1,seq[s+1]])
	u0[2] += residualCapacity * (prob.P[2,seq[s+1]]/prob.W[1,seq[s+1]])
	
	u1 = sol.z + prob.P[:,seq[s]]
	u1[1] -= (prob.W[1,seq[s]] - residualCapacity) * (prob.P[1,seq[s-1]]/prob.W[1,seq[s-1]])
	u1[2] -= (prob.W[1,seq[s]] - residualCapacity) * (prob.P[2,seq[s-1]]/prob.W[1,seq[s-1]])
	
	return u0, u1
end 

function evaluate(weightedObj, y)
	# TODO 
end

# Détermine quel point est obtenu pour la borne supérieure de Martello et Toth 
function chooseBound(prob, weights, u0, u1, iter)
	
	n = size(prob.P)[2]
	
	if domine(u0,u1) 
		push!(upperBound, u0) 
		
	elseif domine(u1,u0) 
		push!(upperBound, u1)
		
	else # Pas de dominance entre u0 et u1
		λeq = (u1[2] - u0[2])/(u0[1] - u0[2] - u1[1] + u1[2]) 
		
		# Le λeq d'égalité est plus petit que le λ critique suivant
		if iter < n && λeq < weights[iter+1] 
			weightedObj = weightedSum(prob, weights[iter+1]) 
			
			# La somme pondérée avec u0 est plus grande lorsque λ > λeq
			if evaluate(weightedObj, u0) >= evaluate(weightedObj, u1) 
				push!(upperBound, u0) 
			else 
				push!(upperBound, u1) 
			end 
			
		# Le λeq d'égalité est plus grand que le λ critique précédent
		elseif iter > 1 && λeq > weights[iter-1] 
			weightedObj = weightedSum(prob, weights[iter-1])
			
			# La somme pondérée avec u0 est plus grande lorsque λ < λeq
			if evaluate(weightedObj, u0) <= evaluate(weightedObj, u1) 
				push!(upperBound, u0) 
			else 
				push!(upperBound, u1) 
			end 
		
		else 
			push!(upperBound, u0)
			push!(upperBound, u1)
		end
	end
	
end

function martelloAndToth(prob::_MOMKP) 

	upperBound = []
	
	# Calcul des ratios 
	r1, r2 = ratios(prob)
	
	# Calcul des poids critiques
	weights, pairs = criticalWeights(prob, r1, r2)
	
	# Tri des ratios dans l'ordre lexicographique décroissant selon (r1,r2)
	seq = sortperm(r1, rev=true) # Séquence d'objets
	pos = sortperm(seq)          # Position des objets dans la séquence
	
	# Construction de la première solution dantzig 
	sol, s, residualCapacity = dantzigSolution(prob, seq) 
	
	# Calcul de la borne de Martello et Toth
	u0, u1 = uMT(prob, seq, sol, residualCapacity, s) 
	
		
end
