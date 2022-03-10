################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Branch-and-bound algorithm                                           #
################################################################################

include("branchAndBoundFunctions.jl")
include("lpRelaxation.jl")

# ----- BRANCH-AND-BOUND ----------------------------------------------------- #
# Recursive branching function 
function branch!(η::Node, 
                 prob::_MOMKP, 
                 L::Vector{Solution}, 
                 branchingVariables::Vector{Int}, 
                 depth::Int)

    # Compute the upper bound set for η 
    UBη = parametricMethod(prob, L, η.init, η.solInit) 

    println("\ndepth = ", depth, " L = ", L)

    # Compare with lower bound set and update status 
    if UBη.points == [[0.,0.]]
        η.pruned = INFEASIBILITY
    elseif length(L) > 1 && isDominated(UBη, L) 
        η.pruned = DOMINANCE 
        plotBoundSets(UBη, L)
    end 
    # INFEASIBILITY, OPTIMALITY ? 
    println("status : ", η.pruned)

    if η.pruned == NOTPRUNED

        # Branching 
        var = branchingVariables[depth] 
        newInit = setVariable(η.init, var)
        println("Variable ", var, " has been set")

        # var is set to 0
        η0 = Node(DualBoundSet{Float64}(), η.solInit, newInit, NOTPRUNED)
        branch!(η0, prob, L, branchingVariables, depth+1) 

        # var is set to 1
        if prob.W[1,var] <= solInit.ω_ 
            solInit1 = Solution(η.solInit.X[1:end], 
                        η.solInit.z + prob.P[:,var],
                        η.solInit.ω_ - prob.W[1,var])  
            solInit1.X[var] = 1 

            η1 = Node(DualBoundSet{Float64}(), solInit1, newInit, NOTPRUNED) 
            branch!(η1, prob, L, branchingVariables, depth+1) 
        end 
    end 
end 

function branchAndBound(prob::_MOMKP) 

    init = initialisation(prob)
    solInit = Solution{Float64}(didactic) 
    rootNode = 
        Node(DualBoundSet{Float64}(), solInit, init, NOTPRUNED)

    # Lower bound set 
    L = Solution[] 

    # Computes the ranks for each variable 
    rank1, rank2 = ranks(init.r1, init.r2)
    # Branching strategy 
    branchingVariables = sumRank(rank1, rank2, INCREASING)
    println("Branching strategy : ", branchingVariables)

    # Recursive branching function 
    branch!(rootNode, prob, L, branchingVariables, 1)

    return L 
end 