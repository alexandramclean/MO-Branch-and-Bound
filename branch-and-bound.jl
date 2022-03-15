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
                 L::Vector{Solution{T}}, 
                 branchingVariables::Vector{Int}, 
                 depth::Int,
                 method::Method) where T<:Real
    
    verbose = false
    
    add!(L, η.solInit)

    # Upper bound and dominance test 
    if η.status == NOTPRUNED
        # Compute the upper bound set for η 
        η.UB = parametricMethod(prob, L, η.init, η.solInit) 

        # Compare with lower bound set and update status 
        if η.UB.points == [[0.,0.]]
            η.status = INFEASIBILITY
        elseif length(L) > 1 && isDominated(η.UB, L) 
            η.status = DOMINANCE 
            #plotBoundSets(η.UB, L)
        end
    end 

    if verbose 
        println("\ndepth = ", depth)
        printBoundSet(L)
        println(η.UB.points)
        println("solInit.X = ", η.solInit.X)
        println("solInit.z = ", η.solInit.z)
        println("solInit.ω_ = ", η.solInit.ω_)
        println("status : ", η.status)
    end 

    # Branching 
    if η.status == NOTPRUNED
        
        if depth <= length(branchingVariables)
            # Set variable  
            var = branchingVariables[depth] 
            newInit = setVariable(η.init, var, method)
            if verbose 
                println("Variable ", var, " has been set")
            end 

            # var is set to 1
            if prob.W[1,var] <= η.solInit.ω_ 
                solInit1 = Solution(η.solInit.X[1:end], 
                            η.solInit.z + prob.P[:,var],
                            η.solInit.ω_ - prob.W[1,var])  
                solInit1.X[var] = 1 

                η1 = Node(DualBoundSet{Float64}(), η, solInit1, newInit, NOTPRUNED) 
                branch!(η1, prob, L, branchingVariables, depth+1, method) 
            else 
                η1 = Node(DualBoundSet{Float64}(), η, η.solInit, newInit, INFEASIBILITY)
                branch!(η1, prob, L, branchingVariables, depth+1, method)
            end

            # var is set to 0
            η0 = Node(DualBoundSet{Float64}(), η, η.solInit, newInit, NOTPRUNED)
            branch!(η0, prob, L, branchingVariables, depth+1, method) 
 
        else 
            η.status = MAXDEPTH
            if verbose 
                println("Max depth has been reached")
            end 
        end 
    end 
end 

function branchAndBound(prob::_MOMKP, method::Method=PARAMETRIC_LP) 

    # Lower bound set and initial solution  
    if method == DICHOTOMIC
        L       = Solution{Rational{Int}}[] 
        solInit = Solution{Rational{Int}}(prob)
    else 
        L       = Solution{Float64}[]
        solInit = Solution{Float64}(prob) 
    end 

    # Root node 
    init     = initialisation(prob, method)
    rootNode = 
        Node(DualBoundSet{Float64}(), nothing, solInit, init, NOTPRUNED)

    # Computes the ranks for each variable 
    rank1, rank2 = ranks(init.r1, init.r2)
    # Branching strategy 
    branchingVariables = sumRank(rank1, rank2, INCREASING)
    println("Branching strategy : ", branchingVariables)

    # Recursive branching function 
    branch!(rootNode, prob, L, branchingVariables, 1, method)

    return L 
end 