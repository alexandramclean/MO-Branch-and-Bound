################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Calcul paramétrique de la relaxation continue                        #
################################################################################

include("dataStructures.jl")
include("functions.jl")
# Précondition : cas uni-dimensionnel

# Construction d'une solution 
# S'il n'y a pas d'objet cassé on considère que s = l'indice du premier objet 
# non-inséré
function buildSolution(prob::_MOMKP, sequence::Vector{Int64})

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
		
# Calcul de la relaxation continue
function relaxationContinue(prob::_MOMKP)

	upperBound = Solution[]

	# Calcul des ratios
	r1, r2 = ratios(prob)
	
	# Calcul des poids critiques
	weights, pairs = criticalWeights(prob, r1, r2)
	
	# Tri des ratios dans l'ordre lexicographique décroissant selon (r1,r2)
	seq = sortperm(r1, rev=true) # Séquence d'objets
	pos = sortperm(seq)          # Position des objets dans la séquence
	
	# Construction de la première solution
	sol, s, residualCapacity = buildSolution(prob, seq)
	push!(upperBound, sol)
	
	# Boucle principale
	for iter in 1:length(weights)

		(i,j) = pairs[iter] 
		k = min(pos[i], pos[j])	
		sol = copy(sol)

		# La place de l'objet cassé est échangée avec un objet dans le sac
		# OU La place d'un objet dans le sac est échangée avec un objet qui
		# n'est pas dans le sac
		if k == s-1 
			# Enlever les objets s-1 et s 
			residualCapacity += prob.W[1,seq[s-1]] 
			sol.z -= prob.P[:,seq[s-1]]
			sol.X[seq[s-1]] = 0 
			
			sol.z -= sol.X[seq[s]] * prob.P[:,seq[s]] 
			sol.X[seq[s]] = 0 
			
			# Insérer l'objet s 
			if prob.W[1,seq[s]] <= residualCapacity # s est inséré en entier
				addItem!(prob, sol, seq[s])
				residualCapacity -= prob.W[1,seq[s]] 
				
				if residualCapacity > 0 # Il reste de la place dans le sac 
				# Une fraction de l'objet à la position s-1 est insérée
					addBreakItem!(prob, sol, residualCapacity, seq[s-1])
				end 
				# La position de l'objet cassé ne change pas 
				
			else # L'objet s reste l'objet cassé 
				addBreakItem!(prob, sol, residualCapacity, seq[s]) 
				# La position de l'objet cassé change
				s = s-1 
			end
		
		# Cas plusieurs λ identiques
		if !(iter < length(weights) && weights[iter] == weights[iter+1])	
			push!(upperBound, sol)
		end 
					
		# La place de l'objet cassé est échangée avec un objet qui n'est pas 
		# dans le sac
		elseif k == s 
			# Enlever l'objet à la position s 
			sol.z -= sol.X[seq[s]] * prob.P[:,seq[s]] 
			sol.X[seq[s]] = 0 
			
			if prob.W[1,seq[s+1]] <= residualCapacity 
			# L'objet à la place s+1 peut être inséré en entier
				addItem!(prob, sol, seq[s+1]) 
				residualCapacity -= prob.W[1,seq[s+1]]  
				
				if residualCapacity > 0 # Il reste de la place dans le sac 
				# Une fraction de l'objet à la position s est insérée
					addBreakItem!(prob, sol, residualCapacity, seq[s]) 
				end
				# La position de l'objet cassé change
				s = s+1 
				
			else # L'objet à la position s+1 devient l'objet cassé 
				addBreakItem!(prob, sol, residualCapacity, seq[s+1])
				# La position de l'objet cassé ne change pas
			end 
		
		# Cas plusieurs λ identiques
		if !(iter < length(weights) && weights[iter] == weights[iter+1])	
			push!(upperBound, sol)
		end 
		
		end
			
		# Mise à jour de la séquence
		tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
		seq[pos[i]] = i ; seq[pos[j]] = j
	end
	
	return upperBound		

end


# Fixer une variable à 1
function setVariable!(weights::Vector{Rational{Int}}, 
					  pairs::Vector{Tuple{Int64,Int64}}, 
					  var::Int64)
	for i in length(weights):-1:1 
		if var in pairs[i] 
			deleteat!(weights, i)
			deleteat!(pairs, i)
		end
	end
end


