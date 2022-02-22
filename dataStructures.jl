################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Structures de données                                                #
################################################################################

# Datastructure of a multi-objective multi-dimensionnal KP with 0/1 variables
struct _MOMKP
    P  :: Matrix{Int} # profit of items for the objectives, k=1..p, j=1..n
    W  :: Matrix{Int} # weight of items for the constraints, i=1..m, j=1..n
    ω  :: Vector{Int} # capacity of knapsacks, i=1..m
end

# Solution qui peut contenir des fractions d'objets
mutable struct Solution
    X::Vector{Rational{Int}} # Liste des éléments insérés dans le sac
    z::Vector{<:Real}        # Valeurs pour les fonctions objectif
end
Solution(n) = Solution(zeros(Rational{Int},n), [0,0])

function copySolution(sol::Solution)
    return Solution(sol.X[1:end],
                    sol.z[1:end])
end

# Structure qui stocke une valeur de λ et les transpositions correspondantes
struct Transposition
	λ::Rational{Int}			  # Poids critique
	pairs::Vector{Tuple{Int,Int}} # Liste des paires correspondantes
end
Transposition(λ) = Transposition(λ,[])

# Data structure of a constraint generated whilst computing the upper bound set
struct Constraint 
	λ::Rational{Int}       # Critical weight
	point::Vector{Float64} # Associated point (u1,u2)
end
# La contrainte associée est λz1 + (1-λ)z2 <= λu1 + (1-λ)u2 où point = (u1,u2)

@enum Optimisation Max Min

#@enum Status DOMINANCE OPTIMALITY INFEASIBILITY
