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
            parametricMethod(prob, init, η.setvar) 

        #@timeit to "Upper bound" Lη, is_dominated = 
        #    parametricLPrelaxation(prob, L, init, η.setvar, interrupt)

        # Compare with lower bound set and update status 
        #if length(η.UB.points) == 0 || η.UB.points == [[0.,0.]]
        #    η.status = INFEASIBILITY

        #elseif length(η.UB.points) == 1 
        #    η.status = DOMINANCE 

        #else
        if length(L) > 1 

            @timeit to "Dominance" is_dominated = isDominated(η.UB, L)

            if is_dominated #&& 
                # La borne sup ne "dépasse" pas d'un côté ou de l'autre 
                # de la borne inf 
                # max z1 in L < max z1 in UB(η)
                #!(L[end].z[1] < η.UB.points[1][1]
                # min z1 in UB(η) < min z1 in L
                #|| η.UB.points[end][1] < L[1].z[1])

                η.status = DOMINANCE 
            #else 
                plotBoundSets(η.UB, L)
            end
        end

        graphic ? plotBoundSets(η.UB, L) : nothing

        for sol in Lη
            @timeit to "Ordered List" add!(L, sol)
        end 
    end

    #solInit = initialSolution(prob, η.setvar)
    @timeit to "Ordered List" add!(L, η.solInit)

    if verbose 
        println("\ndepth = ", depth)
        println("status : ", η.status)
        println("L = ", [sol.z for sol in L])
        #println("solInit.X = ", η.solInit.X)
        #println("solInit.z = ", η.solInit.z)
        #println("solInit.ω_ = ", η.solInit.ω_)
        println(η.setvar)
    end     

    # Branching 
    if η.status == NOTPRUNED
        
        if depth <= length(branchingVariables)
            # Set variable  
            var = branchingVariables[depth] 

            # var is set to 1
            verbose ? println(var, " is set to 1") : nothing
            setvar1 = setVariable(init, η.setvar, var, 1, method)  

            if prob.W[1,var] <= η.solInit.ω_ 
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

                η1 = Node(nothing, setvar1, solInit1, NOTPRUNED) 
            else 

                η1 = Node(nothing, setvar1, η.solInit, INFEASIBILITY)
            end
            branch!(η1, prob, L, init, branchingVariables, depth+1, method, interrupt)

            # var is set to 0
            verbose ? println(var, " is set to 0") : nothing
            setvar0 = setVariable(init, η.setvar, var, 0, method)

            if η.solInit.ω_ == 0
                # There is no residual capacity, no new solutions can be 
                # obtained in this branch 
                η0 = Node(nothing, setvar0, η.solInit, INFEASIBILITY)
            else 
                η0 = Node(nothing, setvar0, η.solInit, NOTPRUNED)
            end 
            branch!(η0, prob, L, init, branchingVariables, depth+1, method, interrupt)
            
        else 
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
    @timeit to "Initialisation" init = initialisation(prob, method)
    @timeit to "Initial setvar" setvar = initialSetvar(prob, init, method)
    rootNode = 
        Node(nothing, setvar, Solution{Float64}(prob), NOTPRUNED)

    # Computes the ranks for each variable 
    rank1, rank2 = ranks(init.r1, init.r2)
    # Branching strategy 
    branchingVariables = sumRank(rank1, rank2, INCREASING)
    println("Branching strategy : ", branchingVariables)

    # Recursive branching function 
    branch!(rootNode, prob, L, init, branchingVariables, 1, method, interrupt)

    return L 
end 