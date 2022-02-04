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
    z::Vector{Rational{Int}} # Valeurs pour les fonctions objectif
end
Solution(n) = Solution(zeros(Rational{Int},n), [0,0])

function copy(sol::Solution)
    return Solution(sol.X[1:end],
                    sol.z[1:end])
end

@enum Optimisation Max Min

#@enum Status DOMINANCE OPTIMALITY INFEASIBILITY


