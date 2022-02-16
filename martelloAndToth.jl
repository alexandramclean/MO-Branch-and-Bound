################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Calcul paramétrique de la borne de Martello et Toth                  #
################################################################################

include("dataStructures.jl")
include("functions.jl")
include("listeOrdonnee.jl")

# Borne de Martello et Toth
function u0(prob::_MOMKP, seq::Vector{Int}, sol::Solution, s::Int, ω_::Int) 
	U0 = sol.z + [0,0] 
	U0[1] += residualCapacity/prob.W[1,seq[s+1]] * prob.P[1,seq[s+1]]
	U0[2] += residualCapacity/prob.W[1,seq[s+1]] * prob.P[2,seq[s+1]]
	return U0
end

function u1(prob::_MOMKP, seq::Vector{Int}, sol::Solution, s::Int, ω_::Int)
	U1 = sol.z + prob.P[:,seq[s]]
	U1[1] -= (prob.W[1,seq[s]] - residualCapacity)/prob.W[1,seq[s-1]] * prob.P[1,seq[s-1]]
	U1[2] -= (prob.W[1,seq[s]] - residualCapacity)/prob.W[1,seq[s-1]] * prob.P[2,seq[s-1]]
	return U1
end 

function uMT(prob::_MOMKP, seq, sol, s, residualCapacity)
	
	U0 = u0(prob, seq, sol, s, residualCapacity)
	U1 = u1(prob, seq, sol, s, residualCapacity)
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
		println("λ = ", λeq)
		
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

	upperBound = Vector{Float64}[]
	
	# Calcul des ratios 
	r1, r2 = ratios(prob)
	
	# Calcul des poids critiques
	weights, pairs = criticalWeights(prob, r1, r2)
	
	# Regroupement des λ identiques
	transpositions = transpositionPreprocessing(weights, pairs)
	
	# Tri des ratios dans l'ordre lexicographique décroissant selon (r1,r2)
	seq = sortperm(1000000*r1 + r2, rev=true) # Séquence d'objets
	pos = sortperm(seq)          # Position des objets dans la séquence
	
	# Construction de la première solution dantzig 
	sol, s, ω_ = dantzigSolution(prob, seq) 
	
	# Calcul de la première borne de Martello et Toth
	U0, U1 = uMT(prob, seq, sol, s, ω_) 
	
	# Borne à conserver 
	upperBound = Vector{Float64}[]
	chooseBound!(upperBound, weights, U0, U1, 0)
	
	println("U0 = ", U0)
	println("U1 = ", U1)
	
	# Boucle principale
	for iter in 1:length(transpositions) 
	
		if !(length(transpositions[iter].pairs) > 1)
		
			(i,j) = transpositions[iter].pairs[1]
			k = min(pos[i], pos[j])

			if k == s-2
				# Echange des objets s-2 et s-1 

				
			elseif k == s-1 
				# Echange des objets s-1 et s
				
			
			elseif k == s 
				# Echange des objets s et s+1
				
			
			elseif k == s+1
				# Echange des objets s+1 et s+2 			
	
	
			end

			# Update the sequence and positions
			tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
			seq[pos[i]] = i ; seq[pos[j]] = j
	
		end
	end
			
	return upperBound
end


