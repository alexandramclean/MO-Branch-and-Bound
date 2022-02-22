################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Calcul paramétrique de la relaxation continue                        #
################################################################################

include("dataStructures.jl")
include("functions.jl")

# The break item is swapped with an item in the knapsack
function swapWithItemInBag(prob::_MOMKP,
						   seq::Vector{Int},
						   sol::Solution,
						   s::Int,
						   residualCapacity::Int)

	# Remove item s-1
	residualCapacity += prob.W[1,seq[s-1]]
	sol.z -= prob.P[:,seq[s-1]]
	sol.X[seq[s-1]] = 0

	# Remove item s
	sol.z -= sol.X[seq[s]] * prob.P[:,seq[s]]
	sol.X[seq[s]] = 0

	if prob.W[1,seq[s]] <= residualCapacity
		# Item s is inserted
		addItem!(prob, sol, seq[s])
		residualCapacity -= prob.W[1,seq[s]]

		if residualCapacity > 0
			# A fraction of item s-1 is inserted 
			addBreakItem!(prob, sol, residualCapacity, seq[s-1])
		end

	else # Item s remains the break item
		addBreakItem!(prob, sol, residualCapacity, seq[s])
		s = s-1 # The position of the break item changes
	end

	return sol, s, residualCapacity

end

# The break item is swapped with an item that is not in the knapsack
function swapWithItemNotInBag(prob::_MOMKP,
						   	  seq::Vector{Int},
						   	  sol::Solution,
						   	  s::Int,
						   	  residualCapacity::Int)

	# The item in position s is removed
	sol.z -= sol.X[seq[s]] * prob.P[:,seq[s]]
	sol.X[seq[s]] = 0

	if prob.W[1,seq[s+1]] <= residualCapacity
		# Item s+1 can be inserted 
		addItem!(prob, sol, seq[s+1])
		residualCapacity -= prob.W[1,seq[s+1]]

		if residualCapacity > 0
			# A fraction of item s is inserted
			addBreakItem!(prob, sol, residualCapacity, seq[s])
		end
		s = s+1 # The position of the break item changes

	else # Item s+1 becomes the break item
		addBreakItem!(prob, sol, residualCapacity, seq[s+1])
	end

	return sol, s, residualCapacity
end

# Fixer une variable à 1 
function setVariable(prob::_MOMKP, 
					 transpositions::Vector{Transposition},
					 seq::Vector{Int},
					 pos::Vector{Int},
					 var::Int)
						
	newTranspositions = Transposition[]
	newSeq = Vector{Int}(undef,length(seq)-1)
	newPos = deepcopy(pos)
	
	# On enlève la variable des transpositions
	for t in transpositions
	
		if length(t.pairs) > 1
		
			swaps = Tuple{Int,Int}[]
			for pair in pairs 
				if !(var in pair) 
					push!(swaps, pair)
				end
			end
			
			if length(swaps) > 0 
				push!(transpositions, Transposition(t.λ, swaps))
			end
		else
		
			if !(var in t.pairs[1])
				push!(newTranspositions, Transposition(t.λ, t.pairs))
			end
		end
		
	end
	
	# On enlève la variable dans la séquence
	inser = 1
	for i in 1:length(seq)
		if seq[i] != var
			newSeq[inser] = seq[i] 
			inser += 1 
		end
	end
	
	# Mise à jour des positions
	# La position des éléments après var dans la séquence diminue de 1
	for p in pos[var]+1:length(pos)
		newPos[seq[p]] = pos[seq[p]] - 1
	end
	
	return newTranspositions, newSeq, newPos
	
end

# LP relaxation
function parametricMethod(prob::_MOMKP,
						  transpositions::Vector{Transposition},
						  seq::Vector{Int},
						  pos::Vector{Int})

	# Builds the initial solution
	sol, s, ω_ = buildSolution(prob, seq)
	upperBound = Vector{Float64}[]
	push!(upperBound, sol.z)

	nbCasEgalite = 0

	for iter in 1:length(transpositions)

		# Multiple identical critical weights λ
		if length(transpositions[iter].pairs) > 1

			nbCasEgalite += 1

			# Positions corresponding to each transposition
			positions = [(min(pos[i], pos[j]), max(pos[i], pos[j])) 
						for (i,j) in transpositions[iter].pairs]
			sort!(positions)
			
			# Identification of the modified subsequences
			start = positions[1][1] ; finish = positions[1][2]
			
			for p in positions[2:end] 
				if p[1] > finish # Start of a new distinct subsequence
					
					# The subsequence is reversed 
					seq[start:finish] = seq[finish:-1:start]
					
					if start <= s && finish >= s # The solution is modified
						sol, s, ω_ = reoptSolution(prob, seq, start, finish, sol, s, ω_)
						if ω_ > 0
							addBreakItem!(prob, sol, ω_, seq[s])
						end
					end
					
					updatePositions!(seq, pos, start, finish)
					
					start = p[1] ; finish = p[2] 
				else 
					finish = p[2] 
				end
			end 
			
			# The subsequence is reversed 
			seq[start:finish] = seq[finish:-1:start]
					
			if start <= s && finish >= s # The solution is modified
				sol, s, ω_ = reoptSolution(prob, seq, start, finish, sol, s, ω_)
				if ω_ > 0
					addBreakItem!(prob, sol, ω_, seq[s])
				end
			end
					
			updatePositions!(seq, pos, start, finish)

			push!(upperBound, sol.z)

		else

			(i,j) = transpositions[iter].pairs[1]
			k = min(pos[i], pos[j])

			if k == s-1

				sol, s, ω_ = swapWithItemInBag(prob, seq, sol, s, ω_)
				push!(upperBound, sol.z)

			elseif k == s

				sol, s, ω_ = swapWithItemNotInBag(prob, seq, sol, s, ω_)
				push!(upperBound, sol.z)

			end

			# Update the sequence and positions
			tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
			seq[pos[i]] = i ; seq[pos[j]] = j

		end
	end

	println("\tNombre de cas d'égalité : ", nbCasEgalite)

	return upperBound

end
