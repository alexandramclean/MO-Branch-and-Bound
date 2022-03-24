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
                 L::PrimalBoundSet{T}, 
                 branchingVariables::Vector{Int}, 
                 depth::Int,
                 method::Method) where T<:Real
    
    verbose = false
    graphic = false 

    # Upper bound and dominance test 
    if η.status == NOTPRUNED
        # Compute the upper bound set for η 
        @timeit to "Upper bound" η.UB = parametricMethod(prob, L, η.init, η.solInit) 

        # Compare with lower bound set and update status 
        if length(η.UB.points) == 0 || η.UB.points == [[0.,0.]]
            η.status = INFEASIBILITY
        elseif length(L.solutions) > 1 

            println("\nL : ")
            afficher(L.solutions)


            shiftedNadirs = localNadirPoints(L.solutions)

            println()
            afficher(L.nadirs)
            afficher(shiftedNadirs)

            @assert L.nadirs ==  shiftedNadirs "The local nadir points are incorrect"

            #@timeit to "Nadirs" L.nadirs = shiftedLocalNadirPoints(L.solutions)

            @timeit to "Dominance" is_dominated = isDominated(η.UB, L.nadirs)

            if !(L.solutions[end].z[1] < η.UB.points[1][1] 
                # max z1 in L < max z1 in UB(η)
                || η.UB.points[end][1] < L.solutions[1].z[1]) && 
                # min z1 in UB(η) < min z1 in L 
                is_dominated

                η.status = DOMINANCE 
                graphic ? plotBoundSets(η.UB, L.solutions) : nothing
            end
        end
    end 
    #=@timeit to "solInit" add!(L, η.solInit)

    println()
    afficher(L.nadirs)
    afficher(shiftedNadirs)=#

    if verbose 
        println("\ndepth = ", depth)
        println("L = ", [sol.z for sol in L.solutions])
        println("UB = ", η.UB.points)
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
            @timeit to "setVar" newInit = setVariable(η.init, var, method)
            verbose ? println("Variable ", var, " has been set") : nothing

            # var is set to 1
            if prob.W[1,var] <= η.solInit.ω_ 
                solInit1 = Solution(η.solInit.X[1:end], 
                            η.solInit.z + prob.P[:,var],
                            η.solInit.ω_ - prob.W[1,var])  
                solInit1.X[var] = 1 

                η1 = Node(DualBoundSet{Float64}(), newInit, solInit1, NOTPRUNED) 
                branch!(η1, prob, L, branchingVariables, depth+1, method) 
            else 
                η1 = Node(DualBoundSet{Float64}(), newInit, η.solInit, INFEASIBILITY)
                branch!(η1, prob, L, branchingVariables, depth+1, method)
            end

            # var is set to 0
            η0 = Node(DualBoundSet{Float64}(), newInit, η.solInit, NOTPRUNED)
            branch!(η0, prob, L, branchingVariables, depth+1, method) 
 
        else 
            η.status = MAXDEPTH
            verbose ? println("Max depth has been reached") : nothing
        end 
    end 
end 

function branchAndBound(prob::_MOMKP, method::Method=PARAMETRIC_LP) 

    # Lower bound set and initial solution  
    if method == DICHOTOMIC
        L       = PrimalBoundSet{Rational{Int}}()
        solInit = Solution{Rational{Int}}(prob)
    else 
        L       = PrimalBoundSet{Float64}()
        solInit = Solution{Float64}(prob) 
    end 

    # Root node 
    init     = initialisation(prob, method)
    rootNode = 
        Node(nothing, init, solInit, NOTPRUNED)

    # Computes the ranks for each variable 
    rank1, rank2 = ranks(init.r1, init.r2)
    # Branching strategy 
    branchingVariables = sumRank(rank1, rank2, INCREASING)
    println("Branching strategy : ", branchingVariables)

    # Recursive branching function 
    branch!(rootNode, prob, L, branchingVariables, 1, method)

    return L 
end 