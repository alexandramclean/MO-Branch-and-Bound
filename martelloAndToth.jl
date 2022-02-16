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
	U0::Vector{Float64} = sol.z + [0,0] 
	U0[1] += ω_/prob.W[1,seq[s+1]] * prob.P[1,seq[s+1]]
	U0[2] += ω_/prob.W[1,seq[s+1]] * prob.P[2,seq[s+1]]
	return U0
end

function u1(prob::_MOMKP, seq::Vector{Int}, sol::Solution, s::Int, ω_::Int)
	U1::Vector{Float64} = sol.z + prob.P[:,seq[s]]
	U1[1] -= (prob.W[1,seq[s]] - ω_)/prob.W[1,seq[s-1]] * prob.P[1,seq[s-1]]
	U1[2] -= (prob.W[1,seq[s]] - ω_)/prob.W[1,seq[s-1]] * prob.P[2,seq[s-1]]
	return U1
end 

function uMT(prob::_MOMKP, seq::Vector{Int}, sol::Solution, s::Int, ω_::Int)
	
	U0 = u0(prob, seq, sol, s, ω_)
	U1 = u1(prob, seq, sol, s, ω_)
	return U0, U1
end 

# Retourne la valeur de la somme pondérée 
function weightedSum(λ::Rational{Int}, y::Vector{Float64})
	return λ*y[1] + (1 - λ)*y[2]
end

# Retourne vrai si la somme pondérée pour λ est plus grande avec U0
function isBigger(U0::Vector{Float64}, 
				  U1::Vector{Float64}, 
				  λ::Rational{Int})
	return weightedSum(λ, U0) >= weightedSum(λ, U1) 
end

# Détermine quel point est obtenu pour la borne supérieure de Martello et Toth 
function chooseBound!(upperBound::Vector{Vector{Float64}}, 
					  weights::Vector{Rational{Int}}, 
					  U0::Vector{Float64}, 
					  U1::Vector{Float64}, 
					  iter::Int)
	
	if domine(U0,U1) 
		ajouter!(upperBound, U0) 
		
	elseif domine(U1,U0) 
		ajouter!(upperBound, U1)
		
	else # Pas de dominance entre U0 et U1
		λeq = (U1[2] - U0[2])/(U0[1] - U0[2] - U1[1] + U1[2]) 
		println("λ = ", λeq)
		
		# Le λeq d'égalité est plus petit ou égale au λ critique suivant
		if iter < length(weights) && λeq <= weights[iter+1] 
		
			# La somme pondérée avec U0 est plus grande lorsque λ > λeq
			if isBigger(U0, U1, weights[iter+1]) 
				ajouter!(upperBound, U0) 
			else 
				ajouter!(upperBound, U1) 
			end 
			
		# Le λeq d'égalité est plus grand ou égal au λ critique précédent
		elseif iter > 1 && λeq >= weights[iter] 
		
			# La somme pondérée avec U0 est plus grande lorsque λ < λeq
			if isBigger(U0, U1, weights[iter]) 
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
	U0::Vector{Float64}, U1::Vector{Float64} = uMT(prob, seq, sol, s, ω_) 
	
	println("U0 = ", U0)
	println("U1 = ", U1)
	
	# Borne à conserver 
	upperBound = Vector{Float64}[]
	chooseBound!(upperBound, weights, U0, U1, 0)
	
	# Boucle principale
	for iter in 1:6 
	
		println("\nIter ", iter)
		println(transpositions[iter])
	
		if !(length(transpositions[iter].pairs) > 1)
		
			(i,j) = transpositions[iter].pairs[1]
			k = min(pos[i], pos[j])

			if k == s-2
				# Echange des objets s-2 et s-1 
				
				# Update the sequence and positions
				tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
				seq[pos[i]] = i ; seq[pos[j]] = j
				
				# Seul U1 est changé
				U1 = u1(prob, seq, sol, s, ω_)
				
				println("U1 = ", U1)
				
				chooseBound!(upperBound, weights, U0, U1, iter)
				
			elseif k == s-1 
				# Echange des objets s-1 et s
				
				# Solution dantzig modifiée 
				# Objet s-1 retiré 
				sol.z -= prob.P[:,seq[s-1]]
				ω_ += prob.W[1,seq[s-1]] 
				sol.X[seq[s-1]] = 0
				
				if prob.W[1,seq[s]] <= ω_ 
					# Objet s inséré 
					addItem!(prob, sol, seq[s])
					ω_ -= prob.W[1,seq[s]] 
				else 
					# La position de l'objet cassé change
					s = s-1
				end 
				
				# Update the sequence and positions
				tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
				seq[pos[i]] = i ; seq[pos[j]] = j
				
				# U0 et U1 modifiés en conséquence
				U0, U1 = uMT(prob, seq, sol, s, ω_)
				
				println("U0 = ", U0)
				println("U1 = ", U1)
				
				chooseBound!(upperBound, weights, U0, U1, iter)	
							
			
			elseif k == s 
				# Echange des objets s et s+1
				
			
			elseif k == s+1
				# Echange des objets s+1 et s+2 
				
				# Update the sequence and positions
				tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
				seq[pos[i]] = i ; seq[pos[j]] = j	
				
				# Seul U0 est modifié 
				U0 = u0(prob, seq, sol, s, ω_)	
				
				println("U0 = ", U0)
				
				chooseBound!(upperBound, weights, U0, U1, iter)	
	
			else 
			
				println("Inchangé")
				# Update the sequence and positions
				tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
				seq[pos[i]] = i ; seq[pos[j]] = j
	
			end
	
		end
	end
	
	println("\n", upperBound)
			
	return upperBound
end


