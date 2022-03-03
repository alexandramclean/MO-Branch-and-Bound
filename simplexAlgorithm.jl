################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Algorithme du simplexe pour le calcul de la relaxation continue      #
################################################################################

include("functions.jl")

# Computes the reduced costs
function V(prob::_MOMKP, 
           i::Int, # Variable 
           j::Int, # Variable 
           k::Int) # Objective function 
    return prob.P[k,j] - prob.P[k,i]*(prob.W[1,j]//prob.W[1,i])
end 

function reducedCosts(prob::_MOMKP, 
                      c::Int) # Basic variable 
    
    n     = size(prob.P)[2] 
    costs = Matrix{Rational{Int}}(undef,2,n)

    for k in 1:2
        for j in 1:n
            costs[k,j] = V(prob, c, j, k)
        end
    end
    return costs
end

# Candidate variables 
function candidateVariables(costs::Matrix{Rational{Int}}, 
                            sol::Solution) 
    
    n = length(sol.X) 
    candidates = Int[] 

    for j in 1:n
        if sol.X[j] == 0 && costs[1,j] < 0 && costs[2,j] > 0 
            push!(candidates, j)
        elseif sol.X[j] == 1 && costs[1,j] > 0 && costs[2,j] < 0 
            push!(candidates, j) 
        end
    end
    return candidates
end

# Initialisaiton 
function simplexInitialisation(prob::_MOMKP)
    r1, r2    = ratios(prob)
    seq       = sortperm(1000000*r1 + r2, rev=true) 
    return seq 
end 

# Simplex algorithm 
function simplex(prob::_MOMKP, seq::Vector{Int})

    #upperBound = Vector{Float64}[]
    upperBound = Solution[]

    # Lexicographically optimal solution for the first objective function 
    sol, s, _ = buildSolution(prob, seq) 
    push!(upperBound, Solution(sol.X, sol.z))

    # The critical objet constitutes an efficient basic variable 
    c = seq[s] 
    
    costs::Matrix{Rational{Int}} = reducedCosts(prob, c)

    candidates = candidateVariables(costs, sol)
    # Stopping criterion : no candidate variables
    stop = (length(candidates) == 0)

    while !stop 

        costRatios::Vector{Rational{Int}} = [costs[2,j]//costs[1,j] for j in candidates]
        perm = sortperm(costRatios)
        j = candidates[perm[1]] 
    
        δ = prob.W[1,c]//prob.W[1,j] 

        if sol.X[j] == 0 

            if δ*sol.X[c] < 1 

                sol.X[j] = δ*sol.X[c] 
                sol.z   += sol.X[j] * prob.P[:,j] 

                sol.z   -= sol.X[c] * prob.P[:,c] 
                sol.X[c] = 0 

                # xj enters the basis in substitution of xc 
                c = j 

            elseif δ*sol.X[c] > 1 

                sol.X[j] = 1 
                sol.z   += prob.P[:,j] 

                sol.X[c] = sol.X[c] - 1//δ 
                sol.z   -= 1//δ * prob.P[:,c] 

                # xc remains in the basis 

            elseif δ*sol.X[c] == 1 

                sol.X[j] = 1
                sol.z   += prob.P[:,j] 

                sol.z   -= sol.X[c] * prob.P[:,c]
                sol.X[c] = 0

                # The solution is integer 
            end 
        else 

            if δ*(1 - sol.X[c]) < 1 

                sol.X[j] = 1 - δ*(1 - sol.X[c]) 
                sol.z   -= δ * (1 - sol.X[c]) * prob.P[:,j] 

                sol.z   -= sol.X[c] * prob.P[:,c]  
                sol.X[c] = 1 
                sol.z   += prob.P[:,c]
                
                # xj enters the basis in substitution of xc 
                c = j 

            elseif δ*(1 - sol.X[c]) > 1 

                sol.X[j] = 0 
                sol.z   -= prob.P[:,j]

                sol.X[c] = sol.X[c] + 1//δ
                sol.z   += 1//δ * prob.P[:,c] 

                # xc remains in the basis 
            else 

                sol.X[j] = 0
                sol.z   -= prob.P[:,j] 

                sol.z -= sol.X[c] * prob.P[:,c] 
                sol.X[c] = 1
                sol.z   += prob.P[:,c] 

                # The solution is integer
            end 
        end 

        push!(upperBound, Solution(sol.X, sol.z))
        
        costs      = reducedCosts(prob, c)
        candidates = candidateVariables(costs, sol)
        stop       = (length(candidates) == 0)
    end 
        
    return upperBound
end