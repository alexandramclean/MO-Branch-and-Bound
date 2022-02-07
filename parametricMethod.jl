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

	upperBound = Solution[]

	# Calcul des ratios et poids critiques
	r1, r2 = ratios(prob)
	weights, pairs = criticalWeights(prob, r1, r2)

	# Tri des ratios dans l'ordre lexicographique décroissant selon (r1,r2)
	seq = sortperm(r1, rev=true) # Séquence d'objets
	pos = sortperm(seq)          # Position des objets dans la séquence

	# Construction de la première solution
	sol, s, residualCapacity = buildSolution(prob, seq)
	push!(upperBound, sol)

	iter = 1
	# Boucle principale
	while iter <= length(weights)

		sol = copy(sol)

		# Cas plusieurs λ identiques
		if iter < length(weights) && weights[iter] == weights[iter+1]

			# On effectue toutes les transpositions associées à λ
			while iter < length(weights) && weights[iter] == weights[iter+1]
				(i,j) = pairs[iter]

				# Mise à jour de la séquence
				tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
				seq[pos[i]] = i ; seq[pos[j]] = j

				iter += 1
			end

			# Construction de la solution associée à la séquence obtenue
			sol, s, residualCapacity = buildSolution(prob, seq)
			push!(upperBound, sol)

		else
			(i,j) = pairs[iter]
			k = min(pos[i], pos[j])

			# La place de l'objet cassé est échangée avec un objet dans le sac
			if k == s-1
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


# Fixer une variable à 1
function setVariable!(weights::Vector{Rational{Int}},
					  pairs::Vector{Tuple{Int,Int}},
					  var::Int)
	for i in length(weights):-1:1
		if var in pairs[i]
			deleteat!(weights, i)
			deleteat!(pairs, i)
		end
	end
end
