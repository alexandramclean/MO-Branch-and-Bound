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

# Updates the positions after reversing the subsequence between deb and fin 
function updatePositions!(seq::Vector{Int}, 
						  pos::Vector{Int}, 
						  deb::Int, 
						  fin::Int) 

	if fin - deb >= 1
		# Swaps the positions of the items at the edge of the subsequence 
		tmp = pos[seq[deb]] 
		pos[seq[deb]] = pos[seq[fin]]
		pos[seq[fin]] = tmp 
		
		# Appel récursif sur la sous-séquence privée des extrémités
		updatePositions!(seq, pos, deb+1, fin-1)
	end
	
end


# LP relaxation
function parametricMethod(prob::_MOMKP)

	# Computes the ratios and critical weights
	r1, r2 = ratios(prob)
	weights, pairs = criticalWeights(prob, r1, r2)

	# Regroupement des λ identiques
	transpositions = transpositionPreprocessing(weights, pairs)

	# The ratios are sorted in decreasing lexicographical order on (r1,r2)
	seq = sortperm(1000000*r1 + r2, rev=true) # Item sequence
	pos = sortperm(seq) # Positions of the items in the sequence

	# Builds the initial solution
	sol, s, ω_ = buildSolution(prob, seq)
	upperBound = Solution[]
	push!(upperBound, sol)

	nbCasEgalite = 0

	for iter in 1:length(transpositions)

		sol = copySolution(sol)

		# Multiple identical critical weights λ
		if length(transpositions[iter].pairs) > 1

			nbCasEgalite += 1
			prev = deepcopy(seq)

			# Positions corresponding to each transposition
			positions = [(min(pos[i], pos[j]), max(pos[i], pos[j])) 
						for (i,j) in transpositions[iter].pairs]
			sort!(positions)
			
			# Identification of the modified subsequences
			deb = positions[1][1] ; fin = positions[1][2]
			
			for p in positions[2:end] 
				if p[1] > fin # Start of a new distinct subsequence
					
					# The subsequence is reversed 
					seq[deb:fin] = seq[fin:-1:deb]
					
					if deb <= s && fin >= s # The solution is modified
						sol, s, ω_ = reoptSolution(prob, prev, seq, deb, sol, s, ω_)
					end
					
					updatePositions!(seq, pos, deb, fin)
					
					deb = p[1] ; fin = p[2] 
				else 
					fin = p[2] 
				end
			end 
			
			# The subsequence is reversed 
			seq[deb:fin] = seq[fin:-1:deb]
					
			if deb <= s && fin >= s # The solution is modified
				sol, s, ω_ = reoptSolution(prob, prev, seq, deb, sol, s, ω_)
			end
					
			updatePositions!(seq, pos, deb, fin)

			push!(upperBound, sol)

		else

			(i,j) = transpositions[iter].pairs[1]
			k = min(pos[i], pos[j])

			if k == s-1

				sol, s, ω_ = swapWithItemInBag(prob, seq, sol, s, ω_)
				push!(upperBound, sol)

			elseif k == s

				sol, s, ω_ = swapWithItemNotInBag(prob, seq, sol, s, ω_)
				push!(upperBound, sol)

			end

			# Update the sequence and positions
			tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
			seq[pos[i]] = i ; seq[pos[j]] = j

		end
	end

	println("\tNombre de cas d'égalité : ", nbCasEgalite)

	return upperBound

end
