################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Dichotomic method to compute efficient supported solutions           #
################################################################################

include("functions.jl")
include("comboJulia.jl")

# Returns the weighted objective defined by parameters λ1 and λ2
function weightedObjective(prob::_MOMKP,
                           λ1::Int64,
                           λ2::Int64)
    return [λ1*prob.P[1,i] + λ2*prob.P[2,i] for i in 1:size(prob.P)[2]]
end

# Transformer une instance _MOMKP en instance pouvant être lue par combo 
function transformInstance(prob::_MOMKP, 
                           obj::Vector{Int64}) # Fonction objectif pondérée 

    n = size(prob.P)[2]
    probCombo = Vector{item}(undef,n)
    for i in 1:n 
        probCombo[i] = item(obj[i], prob.W[1,i], 0, i-1)
    end 
    return probCombo 
end 

# Retourne la solution pour une fonction objectif pondérée donnée 
function getSolution(prob::_MOMKP, obj::Vector{Int64})

    n::Int32 = size(prob.P)[2] 
    ω = prob.ω[1]

    probCombo = transformInstance(prob, obj)
    maxZ = 0 
    for i in 1:n 
        maxZ += probCombo[i].p 
    end 

    val = comboJulia(probCombo, n, ω, maxZ)
    sol = Solution{Rational{Int}}(prob) 

    for i in 1:n 
        index = probCombo[i].i + 1 
        sol.X[index] = probCombo[i].x//1 
        sol.z   += sol.X[index] * prob.P[:,index] 
        sol.ω_  -= sol.X[index] * prob.W[1,index]
    end
    
    return sol
end 

# Calcul des solutions efficaces supportées par méthode dichotomique 
function solveRecursion!(prob::_MOMKP,
                         X_SE::Vector{Solution{Rational{Int}}},
                         x1::Solution{Rational{Int}}, 
                         x2::Solution{Rational{Int}})
    # Calcul de la direction λ
    λ1 = x2.z[2] - x1.z[2]
    λ2 = x1.z[1] - x2.z[1]

    # Fonction objectif pondérée 
    obj::Vector{Int64} = weightedObjective(prob, λ1, λ2)

    # Calcul de la solution
    x = getSolution(prob, obj)
    add!(X_SE, x)

    # Si le point n'est pas sur le segment z(x1)z(x2) on continue la recherche
    if λ1*x.z[1] + λ2*x.z[2] > λ1*x1.z[1] + λ2*x1.z[2]
        solveRecursion!(prob, X_SE, x1, x)
        solveRecursion!(prob, X_SE, x, x2)
    end
end

function dichotomicMethod(prob::_MOMKP)

    X_SE = Vector{Solution{Rational{Int}}}()

    # Calcul des solutions lexicographiquement optimales
    obj12 = weightedObjective(prob, 1, 0)
    x12   = getSolution(prob, obj12)

    obj21 = weightedObjective(prob, 0, 1)
    x21   = getSolution(prob, obj21)
    
    # Appel récursif
    solveRecursion!(prob, X_SE, x12, x21)

    L = Vector{Solution{Float64}}() 
    for sol in X_SE
        push!(L, Solution{Float64}(sol.X, 
            [Float64(sol.z[1]), Float64(sol.z[2])], sol.ω_))
    end 

    return L
end