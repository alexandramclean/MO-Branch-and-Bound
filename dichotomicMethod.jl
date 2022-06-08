################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Méthode dichotomique pour le calcul de la relaxation continue        #
################################################################################

include("functions.jl")

# Builds a solution including the break item for a given sequence
function buildSolutionDicho(prob::_MOMKP, 
							init::Initialisation, 
							setvar::SetVariables,
							seq::Vector{Int})

	n   = size(prob.P)[2]
	sol = Solution{Rational{Int}}(prob)
	for item in setvar.setToOne
		addItem!(prob, sol, item)
	end 
	i = 1

	while i <= n && prob.W[1,seq[i]] <= sol.ω_ && init.r1[seq[i]] >= 0 
		addItem!(prob, sol, seq[i])
		i += 1
	end

	if sol.ω_ > 0 && init.r1[seq[i]] >= 0
		# A fraction of item s is inserted 
		addBreakItem!(prob, sol, seq[i])
	end

	return sol, i
end

# Returns the weighted objective defined by parameters λ1 and λ2
function weightedObjective(prob::_MOMKP,
						   λ1::Union{Int,Rational{Int}},
						   λ2::Union{Int,Rational{Int}})
	return [λ1*prob.P[1,i] + λ2*prob.P[2,i] for i in 1:size(prob.P)[2]]
end

# Returns the solution maximising the weighted objective defined by λ1 and λ2
function solveWeightedSum(prob::_MOMKP,
						  init::Initialisation,
						  #solInit::Solution{Rational{Int}},
						  λ1::Union{Int,Rational{Int}},
						  λ2::Union{Int,Rational{Int}})

	n = size(prob.P)[2]
	obj = weightedObjective(prob, λ1, λ2)
	r_λ = Vector{Rational{Int}}(undef, n)
	for i in 1:n 
		if init.r1[i] == -1 && init.r2[i] == -1 
			r_λ[i] = -1
		else 
			r_λ[i] = obj[i]//prob.W[1,i]
		end 
	end 
	seq  = sortperm(r_λ, rev=true)
	x, s = buildSolutionDicho(prob, init, seq)
	return x, seq[s]
end

# Returns the lexicographically optimal solutions
function lexicographicSolutions!(prob::_MOMKP,
								 UB::DualBoundSet{Rational{Int}},
								 Lη::Vector{Solution{Rational{Int}}},
								 init::Initialisation,
								 setvar::SetVariables)

	println("\nLexicographic solutions \nUB = ", UB.points)
	
	# Lexicographically optimal solution for z^(1,2) 
	seq12  = sortperm(1000000*init.r1 + init.r2, rev=true) 
	x12, s = buildSolutionDicho(prob, init, setvar, seq12)
	println("x12 : ", x12.z)
	updateBoundSets!(UB, Lη, x12, seq12[s])
	println("UB = ", UB.points)
	
	# Lexicographically optimal solution for z^(2,1) 
	seq21  = sortperm(init.r1 + 1000000*init.r2, rev=true) 
	println(seq21)
	x21, s = buildSolutionDicho(prob, init, setvar, seq21) 
	println("x21 : ", x21.z)
	updateBoundSets!(UB, Lη, x21, seq21[s])
	println("UB = ", UB.points)
	println()

	return x12, x21
end

function solveRecursion!(prob::_MOMKP,
						 UB::DualBoundSet{Rational{Int}},
						 Lη::Vector{Solution{Rational{Int}}},
						 init::Initialisation,
						 setvar::SetVariables,
						 x1::Solution{Rational{Int}}, 
						 x2::Solution{Rational{Int}})
	# Calcul de la direction λ
	λ1 = x2.z[2] - x1.z[2]
	λ2 = x1.z[1] - x2.z[1]

	@assert λ1 != 0 || λ2 != 0 

	# Calcul de la solution
	x, breakItem = solveWeightedSum(prob, init, λ1, λ2)
	updateBoundSets!(UB, Lη, λ1//(λ1 + λ2), x, breakItem)

	# Si le point n'est pas sur le segment z(x1)z(x2) on continue la recherche
	if λ1*x.z[1] + λ2*x.z[2] > λ1*x1.z[1] + λ2*x1.z[2]
		solveRecursion!(prob, UB, Lη, init, setvar, x1, x)
		solveRecursion!(prob, UB, Lη, init, setvar, x, x2)
	end
end

function dichotomicMethod(prob::_MOMKP, 		# Bi01KP instance 
						  init::Initialisation, # Utilities
						  setvar::SetVariables) 

	# Upper bound set 
	UB = DualBoundSet{Rational{Int}}()

	# Integer solutions obtained during the computation of the UBS 
	Lη = Vector{Solution{Rational{Int}}}()

	x12, x21 = lexicographicSolutions!(prob, UB, Lη, init, setvar) 	

	solveRecursion!(prob, UB, Lη, init, setvar, x12, x21)

	computeConstraints!(UB)

	return UB, Lη
end
