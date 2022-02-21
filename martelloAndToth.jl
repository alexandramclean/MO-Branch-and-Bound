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

# Retourne la borne pour laquelle la somme pondérée pour λ est plus grande
function returnBiggest(U0::Vector{Float64}, 
				  	   U1::Vector{Float64}, 
				  	   λ::Rational{Int})
	if weightedSum(λ, U0) >= weightedSum(λ, U1)
		return U0 
	else 
		return U1
	end		 
end

# Détermine quel point est obtenu pour la borne supérieure de Martello et Toth 
# sur l'intervalle [prev, next]
function chooseBound!(upperBound::Vector{Vector{Float64}}, 
					  constraints::Vector{Constraint},
					  prev::Rational{Int}, # Previous λ 
					  next::Rational{Int}, # Next λ 
					  U0::Vector{Float64}, 
					  U1::Vector{Float64})
	
	len = length(upperBound)
	
	if domine(U0,U1) 
		push!(upperBound, U0) 
		if length(upperBound) != len
			push!(constraints, Constraint(prev, U0))
		end
		
	elseif domine(U1,U0) 
		push!(upperBound, U1)
		if length(upperBound) != len
			push!(constraints, Constraint(prev, U1))
		end
		
	else # No dominance between U0 and U1
		λeq = (U1[2] - U0[2])/(U0[1] - U0[2] - U1[1] + U1[2]) 

		if λeq < prev && λeq > next
		
			if weightedSum(prev, U0) >= weightedSum(prev, U1)
				# The weighted sum with U0 is bigger when λ in [λeq, prev]
				push!(upperBound, U0)
				
				if length(upperBound) != len
					push!(constraints, Constraint(prev, U0))
				end
				
				# The weighted sum with U1 is bigger when λ in [next, λeq]
				push!(upperBound, U1)
				
				if length(upperBound) != len
					push!(constraints, Constraint(λeq, U1))
				end
				
			else 
				# The weighted sum with U1 is bigger when λ in [λeq, prev]
				push!(upperBound, U1)
				
				if length(upperBound) != len
					push!(constraints, Constraint(prev, U1))
				end
				
				# The weighted sum with U0 is bigger when λ in [next, λeq]
				push!(upperBound, U0)
				
				if length(upperBound) != len
					push!(constraints, Constraint(λeq, U0))
				end
			end 
			
		else 
		
			if λeq <= next 
				U = returnBiggest(U0, U1, next)
				
			elseif λeq >= prev
				U = returnBiggest(U0, U1, prev)
			end
			
			push!(upperBound, U)
			
			if length(upperBound) != len
				push!(constraints, Constraint(prev, U))
			end
		end
	end
	
end

function martelloAndToth(prob::_MOMKP, 
						 transpositions::Vector{Transposition},
						 seq::Vector{Int},
						 pos::Vector{Int}) 

	upperBound  = Vector{Float64}[]
	constraints = Constraint[] 
	
	# Construction de la première solution dantzig 
	sol, s, ω_ = dantzigSolution(prob, seq) 
	
	# Calcul de la première borne de Martello et Toth
	U0, U1 = uMT(prob, seq, sol, s, ω_) 
	
	# Borne à conserver 
	chooseBound!(upperBound, constraints, 1//1, transpositions[1].λ, U0, U1)

	nbCasEgalite = 0
	
	# Boucle principale
	for iter in 1:length(transpositions)
	
		# Poids critiques précédents et suivants
		prev = transpositions[iter].λ
		if iter == length(transpositions)
			next = 0//1
		else 
			next = transpositions[iter+1].λ
		end
	
		# Cas d'égalité
		if length(transpositions[iter].pairs) > 1
			
			nbCasEgalite += 1
			
			# Positions corresponding to each transposition
			positions = [(min(pos[i], pos[j]), max(pos[i], pos[j])) 
						for (i,j) in transpositions[iter].pairs]
			sort!(positions)
			
			# Identification of the modified subsequences
			subsequences = Tuple{Int,Int}[] 
			start = positions[1][1] ; finish = positions[1][2]
			
			for p in positions[2:end] 
				if p[1] > finish # Start of a new distinct subsequence
					push!(subsequences, (start, finish))
					start = p[1] ; finish = p[2] 	
				else 
					finish = p[2] 
				end
			end 
			push!(subsequences, (start, finish))
			
			# Reversing the subsequences
			for (start, finish) in subsequences
			
				# The subsequence is reversed 
				seq[start:finish] = seq[finish:-1:start]	
				updatePositions!(seq, pos, start, finish)
			
				if start < s-1 && finish == s-1 		# Seul U1 est modifié 
					U1 = u1(prob, seq, sol, s, ω_)	
						
				elseif start == s+1 && finish > s+1		# Seul U0 est modifié
					U0 = u0(prob, seq, sol, s, ω_)
						
				elseif start <= s && finish >= s 
					# La solution dantzig est potentiellement modifiée
					sol, s, ω_ = reoptSolution(prob, seq, start, finish, sol, s, ω_)
					U0, U1 = uMT(prob, seq, sol, s, ω_)

				end
				
			end

			chooseBound!(upperBound, constraints, prev, next, U0, U1)	
			
		else
		
			(i,j) = transpositions[iter].pairs[1]
			k = min(pos[i], pos[j])

			if k == s-2
				# Echange des objets s-2 et s-1 
				
				# Update the sequence and positions
				tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
				seq[pos[i]] = i ; seq[pos[j]] = j
				
				# Seul U1 est changé
				U1 = u1(prob, seq, sol, s, ω_)
				
				chooseBound!(upperBound, constraints, prev, next, U0, U1)
				
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
				
				chooseBound!(upperBound, constraints, prev, next, U0, U1)
								
			elseif k == s 
				# Echange des objets s et s+1
				
				# Solution dantzig potentiellement modifiée 
				if prob.W[1,seq[s+1]] <= ω_
					addItem!(prob, sol, seq[s+1])
					ω_ -= prob.W[1,seq[s+1]]
				end
				
				# Update the sequence and positions
				tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
				seq[pos[i]] = i ; seq[pos[j]] = j
				
				# U0 et U1 modifiés
				U0, U1 = uMT(prob, seq, sol, s, ω_)
				
				chooseBound!(upperBound, constraints, prev, next, U0, U1)
			
			elseif k == s+1
				# Swap items s+1 et s+2 
				
				# Update the sequence and positions
				tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
				seq[pos[i]] = i ; seq[pos[j]] = j	
				
				# Only U0 is modified 
				U0 = u0(prob, seq, sol, s, ω_)	
				
				chooseBound!(upperBound, constraints, prev, next, U0, U1)	
				
			else 

				# Update the sequence and positions
				tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
				seq[pos[i]] = i ; seq[pos[j]] = j
	
			end
	
		end
	end
	
	println("\tNombre de cas d'égalité : ", nbCasEgalite)
	return upperBound, constraints
end


