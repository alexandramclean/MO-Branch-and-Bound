################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Méthode dichotomique pour le calcul de la relaxation continue        #
################################################################################

include("functions.jl")

# Builds a solution including the break item for a given sequence
function buildSolutionDicho(prob::_MOMKP, 
							init::Initialisation, 
							seq::Vector{Int})

	n   = size(prob.P)[2]
	ω_  = prob.ω[1]
	sol = Solution{Rational{Int}}(prob)
	i   = 1

	while i <= n && prob.W[1,seq[i]] <= sol.ω_ && init.r1[i] >= 0 
		item = seq[i]
		# L'objet est inséré
		addItem!(prob, sol, item)
		i += 1
	end

	if ω_ > 0 && init.r1[i] >= 0
		# Une fraction de l'objet s est insérée
		addBreakItem!(prob, sol, seq[i])
	end

	return sol, i
end

# Retourne une fonction objectif pondérée définie par les paramètres λ1 et λ2
function weightedObjective(prob::_MOMKP,
						   λ1::Union{Int,Rational{Int}},
						   λ2::Union{Int,Rational{Int}})
	return [λ1*prob.P[1,i] + λ2*prob.P[2,i] for i in 1:size(prob.P)[2]]
end

# Retourne la solution obtenue en maximisant la somme pondérée donnée par λ1 et λ2
function solveWeightedSum(prob::_MOMKP,
						  init::Initialisation,
						  λ1::Union{Int,Rational{Int}},
						  λ2::Union{Int,Rational{Int}})

	n = size(prob.P)[2]
	obj = weightedObjective(prob, λ1, λ2)
	r_λ = Vector{Rational{Int}}(undef, n)
	for i in 1:n 
		if init.r1[i] == -1 && init.r2[i] == -1 
			r_λ[i] = -1
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
								 UB::DualBoundSet,
								 L::Vector{Solution{Rational{Int}}},
								 init::Initialisation)
	
	# Lexicographically optimal solution for (1,2) 
	seq12  = sortperm(1000000*r1 + r2, rev=true) 
	x12, s = buildSolutionDicho(prob, init, seq12)
	updateBoundSets!(UB, L, 1//1, x12, seq12[s])
	
	# Lexicographically optimal solution for (2,1) 
	seq21  = sortperm(r1 + 1000000*r2, rev=true) 
	x21, s = buildSolutionDicho(prob, init, seq21) 
	updateBoundSets!(UB, L, 0//1, x21, seq21[s])

	return x12, x21
end

function solveRecursion!(prob::_MOMKP,
						 UB::DualBoundSet,
						 L::Vector{Solution{T}},
						 init::Initialisation,
						 x1::Solution{Rational{Int}}, 
						 x2::Solution{Rational{Int}}) where T<:Real
	# Calcul de la direction λ
	λ1 = x2.z[2] - x1.z[2]
	λ2 = x1.z[1] - x2.z[1]

	# Calcul de la solution
	x, breakItem = solveWeightedSum(prob, init, λ1, λ2)
	updateBoundSets!(UB, L, λ1//(λ1 + λ2), x, breakItem)

	# Si le point n'est pas sur le segment z(x1)z(x2) on continue la recherche
	if λ1*x.z[1] + λ2*x.z[2] > λ1*x1.z[1] + λ2*x1.z[2]
		solveRecursion!(prob, UB, L, init, x1, x)
		solveRecursion!(prob, UB, L, init, x, x2)
	end
end

function dichotomicMethod(prob::_MOMKP,
						  L::Vector{Solution{T}},
						  init::Initialisation) where T<:Real

	n = size(prob.P)[2]

	# Upper bound set 
	UB = DualBoundSet{Rational{Int}}()

	# Calcul des solutions lexicographiquement optimales
	x12, x21 = lexicographicSolutions!(prob, UB, L, init) 	

	# Appel récursif
	solveRecursion!(prob, UB, L, init, x12, x21)

	return UB
end
