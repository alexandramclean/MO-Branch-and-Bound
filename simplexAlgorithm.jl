################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Algorithme du simplexe pour le calcul de la relaxation continue      #
################################################################################

include("functions.jl")

# Groups together equivalent items 
function groupEquivalentItems(prob::_MOMKP)

    n      = size(prob.P)[2]
    r1, r2 = ratios(prob) 
    P1     = Int[] 
    P2     = Int[] 
    W      = Int[] 
    done   = Int[] # Items that are equivalent to another item that has already 
    # been processed

    for i in 1:n
        
        if !(i in done)

            p1 = prob.P[1,i]
            p2 = prob.P[2,i]
            w  = prob.W[1,i] 

            for j in i+1:n 

                if r1[i] == r1[j] && r2[i] == r2[j] 
                    # Items i and j are equivalent 
                    p1 += prob.P[1,j] 
                    p2 += prob.P[2,j] 
                    w  += prob.W[1,j] 
                    push!(done, j)
                end 
            end 

            push!(P1, p1)
            push!(P2, p2)
            push!(W, w) 
        end 
    end 

    P = Matrix(undef, 2, length(P1))
    for i in 1:length(P1) 
        P[1,i] = P1[i] ; P[2,i] = P2[i]  
    end 

    return _MOMKP(P, reshape(W, 1, length(W)), prob.ω)
end 

# Computes the reduced costs
function V(prob::_MOMKP, 
           i::Int, # Variable 
           j::Int, # Variable 
           k::Int) # Objective function 
    return prob.P[k,j] - prob.P[k,i]*(prob.W[1,j]//prob.W[1,i])
end 

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
function degeneration(costs::Matrix{Rational{Int}}, 
                      sol::Solution, 
                      c::Int) # Basic variable

    n = length(sol.X) 
    candidates = Int[] 

    for j in 1:n 
        if sol.X[c] == 0 && sol.X[j] == 1 && costs[1,j] > 0 && costs[2,j] < 0 
            push!(candidates, j)
        elseif sol.X[c] == 1 && sol.X[j] == 0 && costs[1,j] < 0 && costs[2,j] > 0
            push!(candidates, j)
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
    println("c = ", c)
    
    costs = reducedCosts(prob, c)
    #println(costs)

    # Candidate variables
    is_integer = isInteger(sol)
    if is_integer 
        candidates = degeneration(costs, sol, c)
    else
        candidates = candidateVariables(costs, sol, c)
    end 

    # Stopping criterion : no candidate variables
    stop = (length(candidates) == 0)
    println("Candidats : ", candidates)

    # Main loop 
    while !stop 

        costRatios = [costs[2,j]//costs[1,j] for j in candidates]
        perm = sortperm(costRatios)
        j = candidates[perm[1]] 
        println("\nCandidat : ", j)
    
        δ = prob.W[1,c]//prob.W[1,j] 
        println("δ = ", δ)

        if sol.X[j] == 0 

            println("cas xj = 0")

            if δ*sol.X[c] < 1 

                sol.X[j] = δ*sol.X[c] 
                sol.z   += sol.X[j] * prob.P[:,j] 

                sol.z   -= sol.X[c] * prob.P[:,c] 
                sol.X[c] = 0 

                # xj enters the basis in substitution of xc 
                c = j 
                is_integer = false

            elseif δ*sol.X[c] > 1 

                sol.X[j] = 1 
                sol.z   += prob.P[:,j] 

                sol.X[c] = sol.X[c] - 1//δ 
                sol.z   -= 1//δ * prob.P[:,c] 

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

            println("cas xj = 1") 

            if δ*(1 - sol.X[c]) < 1 

                sol.X[j] = 1 - δ*(1 - sol.X[c]) 
                sol.z   -= δ * (1 - sol.X[c]) * prob.P[:,j] 

                sol.z   -= sol.X[c] * prob.P[:,c]  
                sol.X[c] = 1 
                sol.z   += prob.P[:,c]
                
                # xj enters the basis in substitution of xc 
                c = j 
                is_integer = false 

            elseif δ*(1 - sol.X[c]) > 1 

                sol.X[j] = 0 
                sol.z   -= prob.P[:,j]

                sol.X[c] = sol.X[c] + 1//δ
                sol.z   += 1//δ * prob.P[:,c] 

                # xc remains in the basis 
                is_integer = false 

            else 

                sol.X[j] = 0
                sol.z   -= prob.P[:,j] 

                sol.z -= sol.X[c] * prob.P[:,c] 
                sol.X[c] = 1
                sol.z   += prob.P[:,c] 
                
                is_integer = true # The solution is integer 
            end 
        end 

        push!(upperBound, sol.z)
        

        costs = reducedCosts(prob, c)

        if is_integer 
            candidates = degeneration(costs, sol, c)
        else 
            candidates = candidateVariables(costs, sol, c) 
        end 
        stop = (length(candidates) == 0)

        println("\nX = ", sol.X)
        println("isInteger : ", is_integer)
        println("z = ", sol.z)
        println("c = ", c)
        println("Candidats : ", candidates)
    end 
        
    return upperBound
end