################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Méthode dichotomique pour le calcul de la relaxation continue        #
################################################################################

include("functions.jl")

# Builds a solution including the break item for a given sequence
function buildSolutionD(prob::_MOMKP, seq::Vector{Int})

	n   = size(prob.P)[2]
	ω_  = prob.ω[1]
	sol = SolutionD(zeros(Rational{Int}, n), [0//1, 0//1])
	i   = 1

	while i <= n && prob.W[1,seq[i]] <= ω_
		item = seq[i]
		# L'objet est inséré
		addItem!(prob, sol, item)
		ω_ -= prob.W[1,item]
		i += 1
	end

	if ω_ > 0
		# Une fraction de l'objet s est insérée
		addBreakItem!(prob, sol, ω_, seq[i])
	end

	return sol, i, ω_
end

# Retourne une fonction objectif pondérée définie par les paramètres λ1 et λ2
function weightedObjective(prob::_MOMKP,
						   λ1::Union{Int,Rational{Int}},
						   λ2::Union{Int,Rational{Int}})
	return [λ1*prob.P[1,i] + λ2*prob.P[2,i] for i in 1:size(prob.P)[2]]
end

# Retourne la solution obtenue en maximisant la somme pondérée donnée par λ1 et λ2
function solveWeightedSum(prob::_MOMKP,
						  λ1::Union{Int,Rational{Int}},
						  λ2::Union{Int,Rational{Int}})

	n = size(prob.P)[2]
	obj = weightedObjective(prob, λ1, λ2)
	r_λ = [obj[i]//prob.W[1,i] for i in 1:n]

	seq = sortperm(r_λ, rev=true)
	x, s, _ = buildSolutionD(prob, seq)
	return x, seq[s]
end

# Returns the lexicographically optimal solutions
function lexicographicSolutions!(prob::_MOMKP,
								 UB::DualBoundSet,
								 L::Vector{Solution},
								 r1::Vector{Rational{Int}},
								 r2::Vector{Rational{Int}}) 
	
	# Lexicographically optimal solution for (1,2) 
	seq12 = sortperm(1000000*r1 + r2, rev=true) 
	x12, s, _ = buildSolutionD(prob, seq12)
	updateBoundSets!(UB, L, x12, seq12[s])
	
	# Lexicographically optimal solution for (2,1) 
	seq21 = sortperm(r1 + 1000000*r2, rev=true) 
	x21, s, _ = buildSolutionD(prob, seq21) 
	updateBoundSets!(UB, L, x21, seq21[s])

	return x12, x21
end

function solveRecursion!(prob::_MOMKP,
						 UB::DualBoundSet,
						 L::Vector{Solution},
						 x1::SolutionD, x2::SolutionD)
	# Calcul de la direction λ
	λ1 = x2.z[2] - x1.z[2]
	λ2 = x1.z[1] - x2.z[1]

	# Calcul de la solution
	x, breakItem = solveWeightedSum(prob, λ1, λ2)
	updateBoundSets!(UB, L, x, breakItem)

	# Si le point n'est pas sur le segment z(x1)z(x2) on continue la recherche
	if λ1*x.z[1] + λ2*x.z[2] > λ1*x1.z[1] + λ2*x1.z[2]
		solveRecursion!(prob, UB, L, x1, x)
		solveRecursion!(prob, UB, L, x, x2)
	end
end

function dichotomicMethod(prob::_MOMKP,
						  L::Vector{Solution},
						  r1::Vector{Rational{Int}},
						  r2::Vector{Rational{Int}})

	n = size(prob.P)[2]

	# Upper bound set 
	UB = DualBoundSet()

	# Calcul des solutions lexicographiquement optimales
	x12, x21 = lexicographicSolutions!(prob, UB, L, r1, r2) 	

	# Appel récursif
	solveRecursion!(prob, UB, L, x12, x21)

	return UB
end
