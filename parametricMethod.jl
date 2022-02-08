################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Calcul paramétrique de la relaxation continue                        #
################################################################################

include("dataStructures.jl")
include("functions.jl")
# Précondition : cas uni-dimensionnel
			
# Calcul de la relaxation continue
function parametricMethod(prob::_MOMKP)

	# Calcul des ratios et poids critiques
	r1, r2 = ratios(prob)
	weights, pairs = criticalWeights(prob, r1, r2)
	
	# Regroupement des λ identiques
	transpositions = transpositionPreprocessing(weights, pairs)		

	# Tri des ratios dans l'ordre lexicographique décroissant selon (r1,r2)
	seq = sortperm(1000*r1 + r2, rev=true) # Séquence d'objets
	pos = sortperm(seq)          		   # Position des objets dans la séquence
	
	println(seq, "\n")
	#=for t in transpositions 
		println("λ = ", t.λ, "\t\t", t.pairs) 
	end=# 
	
	# Construction de la première solution
	sol, s, residualCapacity = buildSolution(prob, seq)
	println("s = ", s)
	
	upperBound = Solution[]
	push!(upperBound, sol)

	iter = 1
	# Boucle principale
	while iter <= length(transpositions)

		#println("Iter ", iter)
		
		sol = copySolution(sol)

		# Cas plusieurs λ identiques
		if length(transpositions[iter].pairs) > 1 

			# On effectue toutes les transpositions associées à λ
			nbTransp = 1
			while nbTransp <= length(transpositions[iter].pairs)
				(i,j) = transpositions[iter].pairs[nbTransp]

				# Mise à jour de la séquence
				tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
				seq[pos[i]] = i ; seq[pos[j]] = j

				nbTransp += 1
			end

			# Construction de la solution associée à la séquence obtenue
			sol, s, residualCapacity = buildSolution(prob, seq)
			push!(upperBound, sol)
			iter += 1

		else
		
			(i,j) = transpositions[iter].pairs[1]
			k = min(pos[i], pos[j]) 

			# La place de l'objet cassé est échangée avec un objet dans le sac
			if k == s-1
			
				@assert (pos[i] == pos[j]+1 || pos[j] == pos[i]+1) "La transposition doit être entre deux éléments successifs de la séquence"
					
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
				
			push!(upperBound, sol)

			# L'objet cassé est échangée avec un objet qui n'est pas dans le sac
			elseif k == s

				@assert (pos[i] == pos[j]+1 || pos[j] == pos[i]+1) "La transposition doit être entre deux éléments successifs de la séquence"
				
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

			push!(upperBound, sol)
			
			end 
			
			# Mise à jour de la séquence
			tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
			seq[pos[i]] = i ; seq[pos[j]] = j

			iter += 1
		end
	end

	return upperBound

end

