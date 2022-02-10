################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Calcul paramétrique de la relaxation continue                        #
################################################################################

include("dataStructures.jl")
include("functions.jl")
# Précondition : cas uni-dimensionnel
	
# La place de l'objet cassé est échangée avec un objet dans le sac	
function swapWithItemInBag(prob::_MOMKP, 
						   seq::Vector{Int}, 
						   sol::Solution, 
						   s::Int, 
						   residualCapacity::Union{Int,Rational{Int}})

	# Enlever les objets s-1 et s
	residualCapacity += prob.W[1,seq[s-1]]
	sol.z -= prob.P[:,seq[s-1]]
	sol.X[seq[s-1]] = 0

	sol.z -= sol.X[seq[s]] * prob.P[:,seq[s]]
	sol.X[seq[s]] = 0

	if prob.W[1,seq[s]] <= residualCapacity
		# L'objet s est inséré en entier
		addItem!(prob, sol, seq[s])
		residualCapacity -= prob.W[1,seq[s]]

		if residualCapacity > 0
			# Une fraction de l'objet s-1 est insérée
			addBreakItem!(prob, sol, residualCapacity, seq[s-1])
		end

	else # L'objet s reste l'objet cassé
		addBreakItem!(prob, sol, residualCapacity, seq[s])
		s = s-1 # La position de l'objet cassé change
	end
	
	return sol, s, residualCapacity 
	
end

# L'objet cassé est échangée avec un objet qui n'est pas dans le sac
function swapWithItemNotInBag(prob::_MOMKP, 
						   	  seq::Vector{Int}, 
						   	  sol::Solution, 
						   	  s::Int, 
						   	  residualCapacity::Union{Int,Rational{Int}})
						   	  
	# L'objet à la position s est enlevé
	sol.z -= sol.X[seq[s]] * prob.P[:,seq[s]]
	sol.X[seq[s]] = 0

	if prob.W[1,seq[s+1]] <= residualCapacity
		# L'objet s+1 peut être inséré en entier
		addItem!(prob, sol, seq[s+1])
		residualCapacity -= prob.W[1,seq[s+1]]

		if residualCapacity > 0
			# Une fraction de l'objet s est insérée
			addBreakItem!(prob, sol, residualCapacity, seq[s])
		end
		s = s+1 # La position de l'objet cassé change

	else # L'objet s+1 devient l'objet cassé
		addBreakItem!(prob, sol, residualCapacity, seq[s+1])
	end
	
	return sol, s, residualCapacity
end
			
# Calcul de la relaxation continue
function parametricMethod(prob::_MOMKP)

	# Calcul des ratios et poids critiques
	r1, r2 = ratios(prob)
	weights, pairs = criticalWeights(prob, r1, r2)
	
	# Regroupement des λ identiques
	transpositions = transpositionPreprocessing(weights, pairs)		

	# Tri des ratios dans l'ordre lexicographique décroissant selon (r1,r2)
	seq = sortperm(1000000*r1 + r2, rev=true) # Séquence d'objets
	pos = sortperm(seq)          		      # Position des objets dans la séquence
	
	# Construction de la première solution
	sol, s, residualCapacity = buildSolution(prob, seq)
	
	upperBound = Solution[]
	push!(upperBound, sol)

	# Boucle principale
	for iter in 1:length(transpositions)
		
		sol = copySolution(sol)

		# Cas plusieurs λ identiques
		if length(transpositions[iter].pairs) > 1 
	
			#= Donne les positions correspondant à chaque transposition à faire 
			positions = [(min(pos[i], pos[j]), max(pos[i], pos[j])) for (i,j) in transpositions[iter].pairs] 
	
			# Trie les positions dans l'ordre croissant
			increasing = transpositions[iter].pairs[sortperm(positions)] 
			
			if checkTranspositions(seq, increasing) 
				# On applique les transpositions dans l'ordre croissant 
				swaps = increasing
			else 
				# On applique les transpositions dans l'ordre décroissant 
				swaps = transpositions[iter].pairs[sortperm(positions, rev=true)] 
			end 
				
			for t in 1:length(swaps) 
					
				(i,j) = swaps[t]
				#=k = min(pos[i], pos[j]) 
					
				if k == s-1 
					sol, s, residualCapacity = swapWithItemInBag(prob, seq, sol, s, residualCapacity)
				elseif k == s 
					sol, s, residualCapacity = swapWithItemNotInBag(prob, seq, sol, s, residualCapacity) 
				end =#
				
				# Mise à jour de la séquence
				tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
				seq[pos[i]] = i ; seq[pos[j]] = j
			end =#
			
			# Donne les positions correspondant à chaque transposition à faire 
			positions = [(min(pos[i], pos[j]), max(pos[i], pos[j])) for (i,j) in transpositions[iter].pairs] 
	
			# Trie les positions dans l'ordre croissant
			increasing = transpositions[iter].pairs[sortperm(positions)] 
	
			if checkTranspositions(seq, increasing) 
				swaps = increasing
			else 
				swaps = transpositions[iter].pairs[sortperm(positions, rev=true)]
			end 
	
			for t in 1:length(swaps) 
				(i,j) = swaps[t] 
			
				# Mise à jour de la séquence
				tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
				seq[pos[i]] = i ; seq[pos[j]] = j
			end
			
			sol, s, residualCapacity = buildSolution(prob, seq)			
			push!(upperBound, sol)
				
		else
		
			(i,j) = transpositions[iter].pairs[1]
			k = min(pos[i], pos[j]) 

			# La place de l'objet cassé est échangée avec un objet dans le sac
			if k == s-1
					
				sol, s, residualCapacity = swapWithItemInBag(prob, seq, sol, s, residualCapacity)
				push!(upperBound, sol)

			# L'objet cassé est échangée avec un objet qui n'est pas dans le sac
			elseif k == s
				
				sol, s, residualCapacity = swapWithItemNotInBag(prob, seq, sol, s, residualCapacity)
				push!(upperBound, sol)
			
			end 
			
			# Mise à jour de la séquence
			tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
			seq[pos[i]] = i ; seq[pos[j]] = j

		end
	end

	return upperBound

end

