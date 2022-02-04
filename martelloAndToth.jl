################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Calcul paramétrique de la borne de Martello et Toth                  #
################################################################################

include("dataStructures.jl")
include("functions.jl")
include("listeOrdonnee.jl")

# Borne de Martello et Toth
function u0(prob::_MOMKP, seq, sol, residualCapacity, s) 
	U0 = sol.z + [0,0] 
	U0[1] += residualCapacity * (prob.P[1,seq[s+1]]/prob.W[1,seq[s+1]])
	U0[2] += residualCapacity * (prob.P[2,seq[s+1]]/prob.W[1,seq[s+1]])
	return U0
end

function u1(prob::_MOMKP, seq, sol, residualCapacity, s)
	U1 = sol.z + prob.P[:,seq[s]]
	U1[1] -= (prob.W[1,seq[s]] - residualCapacity) * (prob.P[1,seq[s-1]]/prob.W[1,seq[s-1]])
	U1[2] -= (prob.W[1,seq[s]] - residualCapacity) * (prob.P[2,seq[s-1]]/prob.W[1,seq[s-1]])
	return U1
end 

function uMT(prob::_MOMKP, 
			 seq,
	     	 sol, # Solution dantzig
	     	 residualCapacity, 
	     	 s) # Position de l'objet cassé 
	
	U0 = u0(prob, seq, sol, residualCapacity, s)
	U1 = u1(prob, seq, sol, residualCapacity, s)
	return U0, U1
end 

# Retourne la valeur de la somme pondérée 
function weightedSum(λ, y)
	return λ*y[1] + (1 - λ)*y[2]
end

# Détermine quel point est obtenu pour la borne supérieure de Martello et Toth 
function chooseBound!(upperBound, weights, U0, U1, iter)
	
	if domine(U0,U1) 
		ajouter!(upperBound, U0) 
		
	elseif domine(U1,U0) 
		ajouter!(upperBound, U1)
		
	else # Pas de dominance entre U0 et U1
		λeq = (U1[2] - U0[2])/(U0[1] - U0[2] - U1[1] + U1[2]) 
		#println("λ = ", λeq)
		
		# Le λeq d'égalité est plus petit que le λ critique suivant
		if iter < length(weights) && λeq < weights[iter+1] 
			# La somme pondérée avec U0 est plus grande lorsque λ > λeq
			if weightedSum(weights[iter+1], U0) >= weightedSum(weights[iter+1], U1) 
				ajouter!(upperBound, U0) 
			else 
				ajouter!(upperBound, U1) 
			end 
			
		# Le λeq d'égalité est plus grand que le λ critique précédent
		elseif iter > 1 && λeq > weights[iter] 
			# La somme pondérée avec U0 est plus grande lorsque λ < λeq
			if weightedSum(weights[iter], U0) >= weightedSum(weights[iter], U1) 
				ajouter!(upperBound, U0) 
			else 
				ajouter!(upperBound, U1) 
			end 
		else 
			ajouter!(upperBound, U0)
			ajouter!(upperBound, U1)
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
	
	# Calcul de la première borne de Martello et Toth
	U0, U1 = uMT(prob, seq, sol, residualCapacity, s) 
	chooseBound!(upperBound, weights, U0, U1, 0)
	
	println("U0 = ", U0)
	println("U1 = ", U1)
	
	# Boucle principale
			
	return upperBound
end

# Exemple didactique
prob = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])


