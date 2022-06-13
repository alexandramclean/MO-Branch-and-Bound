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
                 method::Method,
                 interrupt::Bool) where T<:Real 
    
    verbose = false
    graphic = false 

    # Compute the upper bound set for η 
    if interrupt 
        if method == PARAMETRIC_LP 
            @timeit to "Upper bound" Lη, η.status = 
                parametricLPrelaxation(prob, L, η.init, η.setvar, true)
        else 
            @timeit to "Upper bound" Lη, η.status = 
                dichotomicMethod(prob, L, η.init, η.setvar)
        end 

    elseif method == PARAMETRIC_LP 
        @timeit to "Upper bound" η.UB, Lη = 
            parametricMethod(prob, η.init, η.setvar) 

    elseif method == DICHOTOMIC 
        @timeit to "Upper bound" η.UB, Lη = dichotomicMethod(prob, η.init, η.setvar)
    
    else 
        @timeit to "Upper bound" η.UB, Lη = simplex(prob, η.init, η.setvar)
    end 

    # Pruning 
    if !interrupt 
        @timeit to "prune" η.status = prune(η, L, Lη)
    end 

    # Adding any integer solutions found during the computation of the UBS
    # to the incumbent set 
    for sol in Lη
        @timeit to "Ordered List" add!(L, sol)
    end 

    @timeit to "Ordered List" add!(L, η.solInit)

    if verbose 
        println("depth = ", depth)
        #println("status : ", η.status)
        println("L = ", [sol.z for sol in L])
        println("Lη = ", Lη)
        println("U(η) = ", η.UB.points)
        #println("Constraints : ", η.UB.constraints)
        #println("solInit.X = ", η.solInit.X)
        #println("solInit.z = ", η.solInit.z)
        #println("solInit.ω_ = ", η.solInit.ω_)
        println("setvar = ", η.setvar)
    end     

    # Branching 
    if η.status == NOTPRUNED
        
        if depth <= length(branchingVariables)
            # Set variable  
            var = branchingVariables[depth] 

            # var is set to 1
            verbose ? println("\n", var, " is set to 1") : nothing 

            if prob.W[1,var] <= η.solInit.ω_ 
                # var can be assigned to 1 
                @timeit to "setvar" setvar1 = setVariable(η.setvar, var, 1, method) 

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

                @timeit to "setvarInit" init1 = setVariableInit(η.init, var, method)
                η1 = Node(nothing, setvar1, init1, solInit1, NOTPRUNED) 

                if method == DICHOTOMIC || length(init1.seq) > 0
                    branch!(η1, prob, L, branchingVariables, depth+1, 
                    method, interrupt)
                else 
                    verbose ? println("status : INFEASIBILITY") : nothing 
                end 
            else 
                verbose ? println("status : INFEASIBILITY") : nothing 
            end
            
            # var is set to 0
            verbose ? println("\n", var, " is set to 0") : nothing

            @timeit to "setvar" setvar0 = setVariable(η.setvar, var, 0, method)
            @timeit to "setvarInit" init0 = setVariableInit(η.init, var, method)
            η0 = Node(nothing, setvar0, init0, η.solInit, NOTPRUNED)

            if method == DICHOTOMIC || length(init0.seq) > 0 
                branch!(η0, prob, L, branchingVariables, depth+1, 
                method, interrupt)
            else  
                verbose ? println("status : INFEASIBILITY") : nothing 
            end 

        else 
            # There are no more variables to assign 
            η.status = MAXDEPTH
            verbose ? println("Max depth has been reached") : nothing
        end 
    end 
end 

function branchAndBound(prob::_MOMKP, # Bi01KP instance
                        # Initial lower bound set 
                        L::Vector{Solution{T}}=Vector{Solution{Float64}}(), 
                        # Method for computing the upper bound set 
                        method::Method=PARAMETRIC_LP,   
                        # The computation of the upper bound set can be interrupted              
                        interrupt::Bool=false) where T<:Real

    # Initial solution  
    if method == DICHOTOMIC
        solInit = Solution{Rational{Int}}(prob)
    else 
        solInit = Solution{Float64}(prob) 
    end 

    # Root node 
    @timeit to "Initialisation" init   = initialisation(prob, method)
    @timeit to "Initial setvar" setvar = initialSetvar(prob, init, method)

    # Adding the lexicographically optimal solutions if none are provided 
    if length(L) == 0
        x12, x21 = lexicographicSolutions(prob)

        if method == DICHOTOMIC
            add!(L, x12)
            add!(L, x21)
        else 
            add!(L, Solution{Float64}(x12))
            add!(L, Solution{Float64}(x21))
        end 
    end 

    rootNode = Node(nothing, setvar, init, solInit, NOTPRUNED)

    # Computes the ranks for each variable 
    rank1, rank2 = ranks(init.r1, init.r2)
    # Branching strategy 
    branchingVariables = sumRank(rank1, rank2, INCREASING)

    # Recursive branching function 
    branch!(rootNode, prob, L, branchingVariables, 1, method, interrupt)

    return L 
end 