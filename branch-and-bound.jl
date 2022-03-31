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
                 init::Initialisation, 
                 branchingVariables::Vector{Int}, 
                 depth::Int,
                 method::Method,
                 interrupt::Bool) where T<:Real
    
    verbose = false
    graphic = false

    # Upper bound and dominance test 
    if η.status == NOTPRUNED
        # Compute the upper bound set for η 
        @timeit to "Upper bound" η.UB, Lη = 
            parametricMethod(prob, L, init, η.solInit, interrupt) 

        # Compare with lower bound set and update status 
        if length(η.UB.points) == 0 || η.UB.points == [[0.,0.]]
            η.status = INFEASIBILITY
        elseif length(L.solutions) > 1 

            @timeit to "Dominance" is_dominated = isDominated(η.UB, L)

            if is_dominated && 
                # max z1 in L < max z1 in UB(η)
                !(L.solutions[end].z[1] < η.UB.points[1][1] 
                # min z1 in UB(η) < min z1 in L
                || η.UB.points[end][1] < L.solutions[1].z[1])

                η.status = DOMINANCE 
                graphic ? plotBoundSets(η.UB, L.solutions) : nothing
            end
        end

        for sol in Lη.solutions
            @timeit to "Ordered List" add!(L, sol)
        end 
    end
    @timeit to "Ordered List" add!(L, η.solInit)

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
            verbose ? println("Variable ", var, " has been set") : nothing

            # var is set to 1
            setvar1 = setVariable(init, η.setvar, var, 1, method)  
            weight = 0 
            for var in setvar1.setToOne 
                weight += prob.W[1,var] 
            end 
            if weight <= prob.ω[1] 
                η1 = Node(nothing, setvar1, NOTPRUNED)
            else 
                η1 = Node(nothing, setvar1, INFEASIBILITY)
            end 
            branch!(η1, prob, L, init, branchingVariables, depth+1, method, interrupt)

            #=if prob.W[1,var] <= η.solInit.ω_ 
                if method == DICHOTOMIC
                    @timeit to "solInit" solInit1 = Solution{Rational{Int}}(
                                η.solInit.X[1:end], 
                                η.solInit.z + prob.P[:,var],
                                η.solInit.ω_ - prob.W[1,var])

                else 
                    @timeit to "solInit" solInit1 = Solution{Float64}(
                                η.solInit.X[1:end], 
                                η.solInit.z + prob.P[:,var],
                                η.solInit.ω_ - prob.W[1,var])
                end 
                solInit1.X[var] = 1

                η1 = Node(nothing, newInit, solInit1, NOTPRUNED) 
                branch!(η1, prob, L, branchingVariables, depth+1, method, interrupt) 

            else 
                η1 = Node(nothing, η.init, η.solInit, INFEASIBILITY)
                branch!(η1, prob, L, branchingVariables, depth+1, method, interrupt)
            end=#

            # var is set to 0
            η0 = Node(nothing, 
                      setVariable(init, η.setvar, var, 0, method),
                      NOTPRUNED)
            branch!(η0, prob, L, init, branchingVariables, depth+1, method, interrupt)
            #η0 = Node(nothing, newInit, η.solInit, NOTPRUNED)
            #branch!(η0, prob, L, branchingVariables, depth+1, method, interrupt) 
 
        else 
            η.status = MAXDEPTH
            verbose ? println("Max depth has been reached") : nothing
        end 
    end 
end 

function branchAndBound(prob::_MOMKP, # Bi01KP instance
                        # Initial lower bound set 
                        L::PrimalBoundSet{T}=PrimalBoundSet{Float64}(), 
                        # Method for computing the upper bound set 
                        method::Method=PARAMETRIC_LP,   
                        # The computation of the upper bound set can be interrupted              
                        interrupt::Bool=false) where T<:Real       

    # Lower bound set and initial solution  
    if method == DICHOTOMIC
        #L       = PrimalBoundSet{Rational{Int}}()
        solInit = Solution{Rational{Int}}(prob)
    else 
        #L       = PrimalBoundSet{Float64}()
        solInit = Solution{Float64}(prob) 
    end 

    # Root node 
    @timeit to "Initialisation" init = initialisation(prob, method)
    @timeit to "Initial setvar" setvar = initialSetvar(prob, init, method)
    rootNode = 
        Node(nothing, setvar, NOTPRUNED)

    # Computes the ranks for each variable 
    rank1, rank2 = ranks(init.r1, init.r2)
    # Branching strategy 
    branchingVariables = sumRank(rank1, rank2, INCREASING)
    #println("Branching strategy : ", branchingVariables)

    # Recursive branching function 
    branch!(rootNode, prob, L, init, branchingVariables, 1, method, interrupt)

    return L 
end 