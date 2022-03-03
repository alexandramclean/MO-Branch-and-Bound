################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Data structures                                                      #
################################################################################

# Datastructure of a multi-objective multi-dimensionnal KP with 0/1 variables
struct _MOMKP
    P  :: Matrix{Int} # profit of items for the objectives, k=1..p, j=1..n
    W  :: Matrix{Int} # weight of items for the constraints, i=1..m, j=1..n
    ω  :: Vector{Int} # capacity of knapsacks, i=1..m
end

# Data structure of a solution that can contain fractions of an object
mutable struct Solution
    X::Vector{Rational{Int}} # Vector of binary variables 
    z::Vector{Float64}       # Values for the objective functions 
end
Solution(n) = Solution(zeros(Rational{Int},n), [0,0])

# Solution data structure for the dichotomic method (rational values)
mutable struct SolutionD
    X::Vector{Rational{Int}} # Liste des éléments insérés dans le sac
    z::Vector{Rational{Int}} # Valeurs pour les fonctions objectif
end

# Stores a critical weight λ and the corresponding transposition(s)
struct Transposition
	λ::Rational{Int}			  # Critical weight 
	pairs::Vector{Tuple{Int,Int}} # List of pairs of items to swap 
end
Transposition(λ) = Transposition(λ,[])

# Data structure of a constraint generated whilst computing the upper bound set
struct Constraint 
	λ::Rational{Int}       # Critical weight
	point::Vector{Float64} # Associated point (u1,u2)
end
# The associated constraint is λz1 + (1-λ)z2 <= λu1 + (1-λ)u2

# Data structure representing an upper bound set 
mutable struct DualBoundSet 
    points::Vector{Vector{Float64}}    
    constraints::Vector{Constraint}
    integerSols::Vector{Solution} # Integer solutions
end 
DualBoundSet() = DualBoundSet(Solution[], Constraint[], Int[]) 

@enum Optimisation Max Min

#@enum Status DOMINANCE OPTIMALITY INFEASIBILITY
