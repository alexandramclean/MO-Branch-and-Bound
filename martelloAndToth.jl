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

	U0     = sol.z + [0,0]
	U0[1] += ω_/prob.W[1,seq[s+1]] * prob.P[1,seq[s+1]]
	U0[2] += ω_/prob.W[1,seq[s+1]] * prob.P[2,seq[s+1]]
	return U0
end

function u1(prob::_MOMKP, seq::Vector{Int}, sol::Solution, s::Int, ω_::Int)

	U1     = sol.z + prob.P[:,seq[s]]
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

# Returns the point for which the weighted sum with λ is bigger
function returnBiggest(x::Vector{Float64}, y::Vector{Float64}, λ::Rational{Int})
	weightedSum(λ, x) >= weightedSum(λ, y) ? return x : return y
end

# Determines which point is obtained for the Martello and Toth upper bound in
# the interval [next, prev]
function chooseBound!(upperBound::Vector{Vector{Float64}},
					  constraints::Vector{Constraint},
					  prev::Rational{Int}, # Previous critical weight
					  next::Rational{Int}, # Next critical weight
					  U0::Vector{Float64},
					  U1::Vector{Float64})

	# Constraints to be added if the corresponding point was not already in the
	# upper bound
	len              = length(upperBound)
	constraintsToAdd = Constraint[]

	if domine(U0,U1)
		ajouter!(upperBound, U0)
		push!(constraintsToAdd, Constraint(prev, U0))

	elseif domine(U1,U0)
		ajouter!(upperBound, U1)
		push!(constraintsToAdd, Constraint(prev, U1))

	else # No dominance between U0 and U1

		# λ for which the weighted sums with U0 and U1 are equal
		λeq = (U1[2] - U0[2])/(U0[1] - U0[2] - U1[1] + U1[2])

		if λeq < prev && λeq > next

			if weightedSum(prev, U0) >= weightedSum(prev, U1)
				# The weighted sum with U0 is bigger in [λeq, prev]
				ajouter!(upperBound, U0)
				push!(constraintsToAdd, Constraint(prev, U0))

				# The weighted sum with U1 is bigger in [next, λeq]
				ajouter!(upperBound, U1)
				push!(constraintsToAdd, Constraint(λeq, U1))

			else
				# The weighted sum with U1 is bigger in [λeq, prev]
				ajouter!(upperBound, U1)
				push!(constraintsToAdd, Constraint(prev, U1))

				# The weighted sum with U0 is bigger in [next, λeq]
				ajouter!(upperBound, U0)
				push!(constraintsToAdd, Constraint(λeq, U0))

			end
		else

			if λeq <= next
				U = returnBiggest(U0, U1, next)

			elseif λeq >= prev
				U = returnBiggest(U0, U1, prev)
			end

			ajouter!(upperBound, U)
			push!(constraintsToAdd, Constraint(prev, U))
		end
	end

	if length(upperBound) != len
		for c in constraintsToAdd
			push!(constraints, c)
		end
	end
end

function martelloAndToth(prob::_MOMKP)

	upperBound  = Vector{Float64}[]
	constraints = Constraint[]

	# Calcul des ratios
	r1, r2 = ratios(prob)

	# Calcul des poids critiques
	@time weights, pairs = criticalWeights(prob, r1, r2)

	# Regroupement des λ identiques
	transpositions = transpositionPreprocessing(weights, pairs)

	# Tri des ratios dans l'ordre lexicographique décroissant selon (r1,r2)
	seq = sortperm(1000000*r1 + r2, rev=true) # Item sequence
	pos = sortperm(seq)          			  # Item positions

	# Builds the initial dantzig solution
	sol, s, ω_ = dantzigSolution(prob, seq)

	U0, U1 = uMT(prob, seq, sol, s, ω_)
	chooseBound!(upperBound, constraints, 1//1, weights[1], U0, U1)

	numberCasesIdenticalWeights = 0

	for iter in 1:length(transpositions)

		# Previous and next critical weights
		prev = transpositions[iter].λ
		if iter == length(transpositions)
			next = 0//1
		else
			next = transpositions[iter+1].λ
		end

		# Multiple identical critical weights
		if length(transpositions[iter].pairs) > 1

			numberCasesIdenticalWeights += 1

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

				seq[start:finish] = seq[finish:-1:start]
				updatePositions!(seq, pos, start, finish)

				if start < s-1 && finish == s-1 		# Only U1 is modified
					U1 = u1(prob, seq, sol, s, ω_)

				elseif start == s+1 && finish > s+1		# Only U0 is modified
					U0 = u0(prob, seq, sol, s, ω_)

				elseif start <= s && finish >= s

					# The dantzig solution is potentially modified
					sol, s, ω_ = reoptSolution(prob, seq, start, finish, sol, s, ω_)
					U0, U1 = uMT(prob, seq, sol, s, ω_)
				end
			end

			chooseBound!(upperBound, constraints, prev, next, U0, U1)
		else

			(i,j) = transpositions[iter].pairs[1]
			k = min(pos[i], pos[j])

			# Update the sequence and positions
			tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
			seq[pos[i]] = i ; seq[pos[j]] = j

			if k == s-2 # Swap items s-2 and s-1

				# Only U1 is modified
				U1 = u1(prob, seq, sol, s, ω_)
				chooseBound!(upperBound, constraints, prev, next, U0, U1)

			elseif k == s-1 # Swap items s-1 and s

				# The item previously in position s-1 is removed
				sol.z -= prob.P[:,seq[s]]
				ω_ += prob.W[1,seq[s]]
				sol.X[seq[s]] = 0

				if prob.W[1,seq[s-1]] <= ω_
					# The item previously in position s is inserted
					addItem!(prob, sol, seq[s-1])
					ω_ -= prob.W[1,seq[s-1]]
				else
					# The position of the break item changes
					s = s-1
				end

				U0, U1 = uMT(prob, seq, sol, s, ω_)
				chooseBound!(upperBound, constraints, prev, next, U0, U1)

			elseif k == s # Swap items s and s+1

				if prob.W[1,seq[s]] <= ω_
					# The item previously in position s+1 is inserted
					addItem!(prob, sol, seq[s])
					ω_ -= prob.W[1,seq[s]]
					# The position of the break item changes
					s = s+1
				end

				U0, U1 = uMT(prob, seq, sol, s, ω_)
				chooseBound!(upperBound, constraints, prev, next, U0, U1)

			elseif k == s+1 # Swap items s+1 and s+2

				# Only U0 is modified
				U0 = u0(prob, seq, sol, s, ω_)
				chooseBound!(upperBound, constraints, prev, next, U0, U1)
			end
		end
	end

	println("\tNumber of cases of identical critical weights : ", numberCasesIdenticalWeights)
	return upperBound, constraints
end
