################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Algorithme du simplexe pour le calcul de la relaxation continue      #
################################################################################

include("functions.jl")


function V(prob::_MOMKP, 
           i::Int, # Variable 
           j::Int, # Variable 
           k::Int) # Objective function 
    return prob.P[k,j] - prob.P[k,i]*(prob.W[1,j]//prob.W[1,i])
end 
# Computes the reduced costs
function reducedCosts(prob::_MOMKP, 
                      c::Int) # Basic variable 
    
    n = size(prob.P)[2] 
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
                            sol::Solution, 
                            c::Int) # Basic variable 
    
    n = size(costs)[2]
    candidates = Int[] 

    for j in 1:n
        if j != c 
            if sol.X[j] == 0 && costs[1,j] < 0 && costs[2,j] > 0 
                push!(candidates, j)
            elseif sol.X[j] == 1 && costs[1,j] > 0 && costs[2,j] < 0 
                push!(candidates, j) 
            end 
        end
    end
    return candidates
end

# Degeneration case : the obtained solution is integer
function degeneration(prob::_MOMKP, sol::Solution, c::Int)

    n = size(prob.P)[2] 
    candidates = [] 

    for i in 1:n 
        for j in i+1:n 

            if (sol.X[i] == 0 && sol.X[j] == 1 
                && V(prob, i, j, 1) > 0 && V(prob, i, j, 2) < 0) 
                push!(candidates, (i,j))
            
            elseif (sol.X[i] == 1 && sol.X[j] == 0 
                && V(prob, i, j, 1) < 0 && V(prob, i, j, 2) > 0)
                push!(candidates, (i,j))
            end
        end 
    end 
    return candidates 
end

# Returns true if the solution is integer 
function isInteger(sol::Solution)

    is_integer = true 
    j = 1
    while is_integer && j <= length(sol.X)

        if sol.X[j] > 0 && sol.X[j] < 1 
            is_integer = false 
        end
        j += 1 
    end 
    return is_integer 
end 

# Simplex algorithm 
function simplex(prob::_MOMKP)

    upperBound = Vector{Float64}[]

    # Lexicographically optimal solution for the first objective function 
    r1, r2    = ratios(prob)
    seq       = sortperm(1000000*r1 + r2, rev=true) 
    sol, s, _ = buildSolution(prob, seq) 
    push!(upperBound, sol.z)

    # The critical objet constitutes an efficient basic variable 
    c = seq[s] 
    
    # Compute the reduced costs
    costs = reducedCosts(prob, c)

    # Candidate variables
    is_integer = isInteger(sol)
    if is_integer 
        candidates = degeneration(prob, sol, c)
    else
        candidates = candidateVariables(costs, sol, c)
    end 

    # Stopping criterion : no candidate variables
    stop = (length(candidates) == 0)
    
    # Main loop 
    while !stop 

        costRatios = [costs[2,j]//costs[1,j] for j in candidates]
        j = candidates[sortperm(costRatios)[1]] 
    
        δ = prob.W[1,c]//prob.W[1,j] 
        
        if sol.X[j] == 0 

            if δ*sol.X[c] < 1 

                sol.X[j] = δ*sol.X[c] 
                sol.z   += sol.X[j] * prob.P[:,j] 

                sol.z   -= sol.X[c] * prob.P[:,j] 
                sol.X[c] = 0 

                # xj enters the basis in substitution of xc 
                c = j 
                is_integer = false

            elseif δ*sol.X[c] > 1 

                sol.X[j] = 1 
                sol.z   += prob.P[:,j] 

                sol.X[c] = sol.X[c] - δ 
                sol.z   -= δ * prob.P[:,c] 

                # xc remains in the basis 
                is_integer = false

            elseif δ*sol.X[c] == 1 

                sol.X[j] = 1
                sol.z   += prob.P[:,j] 

                sol.z   -= sol.X[c] * prob.P[:,c]
                sol.X[c] = 0

                is_integer = true # The solution is integer 
            end 
        else 

            if δ*(1 - sol.X[c]) < 1 

                sol.X[j] = 1 - δ*(1 - sol.X[c]) 
                sol.z   -= δ * (1 - sol.X[c]) * prob.P[:,j] 

                sol.X[c] = 1 
                sol.z   += prob.P[:,c]

                c = j # xj enters the basis in substitution of xc 
                is_integer = false 

            elseif δ*(1 - sol.X[c]) > 1 

                sol.X[j] = 0 
                sol.z   -= prob.P[:,j]

                sol.X[c] = sol.X[c] + δ
                sol.z   += δ * prob.P[:,c] 

                # xc remains in the basis 
                is_integer = false 

            else 

                sol.X[j] = 0
                sol.z   -= prob.P[:,j] 

                sol.X[c] = 1
                sol.z   += prob.P[:,c] 
                
                is_integer = true # The solution is integer 
            end 
        end 

        push!(upperBound, sol.z)
        
        costs = reducedCosts(prob, c)
        if is_integer 
            candidates = degeneration(prob, sol, c)
        else 
            candidates = candidateVariables(costs, sol, c) 
        end 

        stop = (length(candidates) == 0)
    end 
        
    return upperBound
end