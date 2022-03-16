################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Calcul paramétrique de la borne de Martello et Toth                  #
################################################################################

include("parametricMethodFunctions.jl")

# Martello and Toth upper bound 
function u0(prob::_MOMKP, seq::Vector{Int}, sol::Solution, s::Int)

	U0::Vector{Float64} = sol.z + [0.0, 0.0]
	U0[1] += sol.ω_/prob.W[1,seq[s+1]] * prob.P[1,seq[s+1]]
	U0[2] += sol.ω_/prob.W[1,seq[s+1]] * prob.P[2,seq[s+1]]
	return U0
end

function u1(prob::_MOMKP, seq::Vector{Int}, sol::Solution, s::Int)

	U1::Vector{Float64} = sol.z + prob.P[:,seq[s]]
	U1[1] -= (prob.W[1,seq[s]] - sol.ω_)/prob.W[1,seq[s-1]] * prob.P[1,seq[s-1]]
	U1[2] -= (prob.W[1,seq[s]] - sol.ω_)/prob.W[1,seq[s-1]] * prob.P[2,seq[s-1]]
	return U1
end

function uMT(prob::_MOMKP, seq::Vector{Int}, sol::Solution, s::Int)

	U0 = u0(prob, seq, sol, s)
	U1 = u1(prob, seq, sol, s)
	return U0, U1
end

# Returns the value of the weighted sum λy1 + (1-λ)y2
function weightedSum(λ::Rational{Int}, y::Vector{Float64})
	return λ*y[1] + (1 - λ)*y[2]
end

# Returns the point for which the weighted sum with λ is bigger
function returnBiggest(x::Vector{Float64}, y::Vector{Float64}, λ::Rational{Int})
	if weightedSum(λ, x) >= weightedSum(λ, y)
		return x 
	else 
		return y
	end
end

# Determines which point is obtained for the Martello and Toth upper bound 
# for λ in the interval [next, prev]
function chooseBound!(UB::DualBoundSet,
					  prev::Rational{Int}, # Previous critical weight
					  next::Rational{Int}, # Next critical weight
					  U0::Vector{Float64},
					  U1::Vector{Float64})

	if dominates(U0,U1)
		updateBoundSet!(UB, prev, U0)

	elseif dominates(U1,U0)
		updateBoundSet!(UB, prev, U1)

	else # No dominance between U0 and U1

		# λ for which the weighted sums with U0 and U1 are equal
		λeq = (U1[2] - U0[2])/(U0[1] - U0[2] - U1[1] + U1[2])

		if λeq < prev && λeq > next

			if weightedSum(prev, U0) >= weightedSum(prev, U1)
				# The weighted sum with U0 is bigger in [λeq, prev]
				updateBoundSet!(UB, prev, U0)

				# The weighted sum with U1 is bigger in [next, λeq]
				updateBoundSet!(UB, λeq, U1)
			else
				# The weighted sum with U1 is bigger in [λeq, prev]
				updateBoundSet!(UB, prev, U1)

				# The weighted sum with U0 is bigger in [next, λeq]
				updateBoundSet!(UB, λeq, U0)
			end
		else

			if λeq <= next
				U = returnBiggest(U0, U1, next)

			elseif λeq >= prev
				U = returnBiggest(U0, U1, prev)
			end

			updateBoundSet!(UB, prev, U)
		end
	end
end

# Computes the Martello and Toth upper bound using the parametric method 
function martelloAndToth(prob::_MOMKP,
						 init::Initialisation,
						 solInit::Solution)

	# Creates copies of the sequence and positions as they will be modified 
	seq = init.seq[1:end] 
	pos = init.pos[1:end] 

	UB  = DualBoundSet{Float64}()

	# Builds the initial dantzig solution
	sol, s = dantzigSolution(prob, seq, solInit)

	U0, U1 = uMT(prob, seq, sol, s)
	chooseBound!(UB, 1//1, init.transpositions[1].λ, U0, U1)

	numberCasesIdenticalWeights = 0

	for iter in 1:length(init.transpositions)

		# Previous and next critical weights
		prev = init.transpositions[iter].λ
		if iter == length(init.transpositions)
			next = 0//1
		else
			next = init.transpositions[iter+1].λ
		end

		# Multiple identical critical weights
		if length(init.transpositions[iter].pairs) > 1

			numberCasesIdenticalWeights += 1

			# Positions corresponding to each transposition
			positions = [(min(pos[i], pos[j]), max(pos[i], pos[j]))
						for (i,j) in init.transpositions[iter].pairs]
			sort!(positions)

			# Identification of the modified subsequences
			subsequences = identifySubsequences(positions)

			# Reversing the subsequences
			for (start, finish) in subsequences

				seq[start:finish] = seq[finish:-1:start]
				updatePositions!(seq, pos, start, finish)

				if start < s-1 && finish == s-1 		# Only U1 is modified
					U1 = u1(prob, seq, sol, s)

				elseif start == s+1 && finish > s+1		# Only U0 is modified
					U0 = u0(prob, seq, sol, s)

				elseif start <= s && finish >= s

					# The dantzig solution is potentially modified
					sol, s = reoptSolution(prob, seq, start, finish, sol)
					U0, U1 = uMT(prob, seq, sol, s)
				end
			end

			chooseBound!(UB, prev, next, U0, U1)
		else

			(i,j) = init.transpositions[iter].pairs[1]
			k = min(pos[i], pos[j])

			# Update the sequence and positions
			tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
			seq[pos[i]] = i ; seq[pos[j]] = j

			if k == s-2     # Swap items s-2 and s-1

				# Only U1 is modified
				U1 = u1(prob, seq, sol, s)
				chooseBound!(UB, prev, next, U0, U1)

			elseif k == s-1 # Swap items s-1 and s

				# The item previously in position s-1 is removed
				sol.z  -= prob.P[:,seq[s]]
				sol.ω_ += prob.W[1,seq[s]]
				sol.X[seq[s]] = 0

				if prob.W[1,seq[s-1]] <= sol.ω_
					# The item previously in position s is inserted
					addItem!(prob, sol, seq[s-1])
				else
					# The position of the break item changes
					s = s-1
				end

				U0, U1 = uMT(prob, seq, sol, s)
				chooseBound!(UB, prev, next, U0, U1)

			elseif k == s   # Swap items s and s+1

				if prob.W[1,seq[s]] <= sol.ω_
					# The item previously in position s+1 is inserted
					addItem!(prob, sol, seq[s])
					# The position of the break item changes
					s = s+1
				end

				U0, U1 = uMT(prob, seq, sol, s)
				chooseBound!(UB, prev, next, U0, U1)

			elseif k == s+1 # Swap items s+1 and s+2

				# Only U0 is modified
				U0 = u0(prob, seq, sol, s)
				chooseBound!(UB, prev, next, U0, U1)
			end
		end
	end

	#println("\tNumber of cases of identical critical weights : ", numberCasesIdenticalWeights)
	return UB
end
