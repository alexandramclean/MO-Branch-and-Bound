################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Algorithme du simplexe pour le calcul de la relaxation continue      #
################################################################################

include("functions.jl")

# Computes the reduced costs 
function reducedCosts(prob::_MOMKP, basicVariable::Int)

end

function simplex(prob::_MOMKP)

    upperBound = Vector{Float64}[]

    # Lexicographically optimal solution for the first objective function 
    r1, r2     = ratios(prob)
    seq        = sortperm(1000000*r1 + r2, rev=true) 
    sol, s, ω_ = buildSolution(prob, seq) 
    push!(upperBound, sol.z)

    # The critical objet constitutes an efficient basic variable 
    basicVariable = seq[s] 
    
    # Compute the reduced costs
    costs = reducedCosts(prob, basicVariable)
    
    # Candidate variables

    # Stopping criterion : no candidate variables

    # Main loop : 
    # Pivot 

    # Degeneration case

    # Reduced costs 

    # Candidate variables

        
end