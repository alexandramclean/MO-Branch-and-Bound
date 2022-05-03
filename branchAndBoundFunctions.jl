################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Branch-and-bound auxiliary functions                                 #
################################################################################

include("functions.jl")
include("displayGraphic.jl")

@enum Order INCREASING DECREASING 

# ----- BRANCHING STRATEGIES ------------------------------------------------- #
# Returns the ranks of items ordered by decreasing order of utility  
function ranks(r1::Vector{Rational{Int}}, # Utilities for objective function 1
               r2::Vector{Rational{Int}}) # Utilities for objective function 2 

    # Ranks for the first objective function 
    rank1 = sortperm(sortperm(r1, rev=true))

    # Ranks for the second objective function 
    rank2 = sortperm(sortperm(r2, rev=true))

    return rank1, rank2 
end 

# Items selected in increasing or decreasing order of min(r1,r2)
function minUtility(r1::Vector{Rational{Int}}, 
                    r2::Vector{Rational{Int}},
                    order::Order)

    min_utility = [min(r1[j], r2[j]) for j in 1:length(r1)]
    if order == DECREASING
        return sortperm(min_utility, rev=true) # Best variable first
    else 
        return sortperm(min_utility)           # Worst variable first
    end
end 

# Items selected in increasing or decreasing order of max(u1,u2) 
function maxUtility(r1::Vector{Rational{Int}}, 
                    r2::Vector{Rational{Int}},
                    order::Order)

    max_utility = [max(r1[j], r2[j]) for j in 1:length(r1)]
    if order == DECREASING 
        return sortperm(max_utility, rev=true) # Best variable first
    else 
        return sortperm(max_utility)           # Worst variable first

    end
end 

# Items selected in increasing or decreasing order of (u1+u2)/2 
function avgUtility(r1::Vector{Rational{Int}}, 
                    r2::Vector{Rational{Int}},
                    order::Order)

    avg_utility = [(r1[j] + r2[j])//2 for j in 1:length(u1)] 
    if order == DECREASING 
        return sortperm(avg_utility, rev=true) # Best variable first
    else 
        return sortperm(avg_utility)           # Worst variable first
    end 
end 

# Items selected in increasing or decreasing order of r1+r2 
function sumRank(rank1::Vector{Int}, 
                 rank2::Vector{Int},
                 order::Order) 
    
    sum_rank = [rank1[j] + rank2[j] for j in 1:length(rank1)]
    if order == INCREASING 
        return sortperm(sum_rank)           # Best variable first
    else 
        return sortperm(sum_rank, rev=true) # Worst variable first
    end 
end 

# Items selected in increasing or decreasing order of min(r1,r2)
function minRank(rank1::Vector{Int},
                 rank2::Vector{Int},
                 order::Order)

    n = length(rank1)
    min_rank = [min(rank1[j],rank2[j]) + (rank1[j] + rank2[j])//2*n for j in 1:n]
    if order == INCREASING 
        return sortperm(min_rank)           # Best variable first
    else
        return sortperm(min_rank, rev=true) # Worst variable first
    end 
end 

# Items selected in increasing or decreasing order of max(r1,r2)
function maxRank(rank1::Vector{Int},
                 rank2::Vector{Int},
                 order::Order)

    n = length(rank1)
    max_rank = [max(rank1[j],rank2[j]) + (rank1[j] + rank2[j])//2*n for j in 1:n]
    if order == INCREASING 
        return sortperm(max_rank)           # Best variable first
    else
        return sortperm(max_rank, rev=true) # Worst variable first
    end 
end 

# ----- LOCAL NADIR POINTS AND DOMINANCE TESTS ------------------------------- # 
# Returns true if the point y verifies all the constraints 
function verifiesConstraints(constraints::Vector{Constraint}, 
                             yl::Vector{Float64},
                             yr::Vector{Float64})
            
    verif = true 
    i = 1 

    # If a shifted local nadir point does not verify one of the constraints
    # the algorithm stops 
    while verif && i <= length(constraints)
        λ     = constraints[i].λ 
        point = constraints[i].point

        if λ!= 1//1
            verif = verif && 
                λ*(yl[1]+1.) + (1-λ)*(yr[2]+1.) <= λ*point[1] + (1-λ)*point[2]
        end 
        i += 1 
    end 
    return verif 
end 

# Returns true if the upper bound set for a particular node is dominated by the 
# lower bound set by using the constraints and shifted local nadir points 
function isDominated(constraints::Vector{Constraint}, 
                     L::Vector{Solution{T}}
                    ) where T<:Real

    is_not_dominated  = false
    i = 2 

    while !is_not_dominated && i <= length(L)
        # If there is a shifter local nadir point that verifies the constraints
        # the node cannot be pruned 
        is_not_dominated = is_not_dominated || 
            verifiesConstraints(constraints, L[i-1].z, L[i].z)
        i += 1 
    end 
    return !is_not_dominated
end 

# ----- PRUNING -------------------------------------------------------------- #
# Determines whether a node is pruned 
# -- Note : The case in which the subproblem is infeasible is managed in the 
# branch-and-bound algorithm 
function prune(η::Node,                # Node 
               L::Vector{Solution{T}}, # Lower bound set 
               # Integer solutions obtained while commputing the upper bound set
               Lη::Vector{Solution{T}}) where T<:Real  

    UB = η.UB.constraints 

    # Optimality : The upper bound is a single, feasible point 
    if length(UB) <= 2 && length(Lη) == 1 && UB[1].point == Lη[1].z 
        return OPTIMALITY
    end 

    # Dominance 
    isDominated = true 

    if length(L) > 1 
        for i in 2:length(L)

            # Local nadir point 
            nadir = [L[i-1].z[1], L[i].z[2]] 
            
            iter = 1
            nadirInUpper = true 

            while iter <= length(UB) && nadirInUpper
                lhs = UB[iter].λ * nadir[1] + (1 - UB[iter].λ) * nadir[2]
                rhs = UB[iter].λ * UB[iter].point[1] + 
                        (1 - UB[iter].λ) * UB[iter].point[2]
                nadirInUpper = nadirInUpper && lhs <= rhs 
                iter += 1 
            end 

            isDominated = isDominated && !nadirInUpper 
        end 

        if isDominated 
            return DOMINANCE 
        else 
            return NOTPRUNED
        end
    else 
        return NOTPRUNED
    end 
end 

# ----- GRAPHIC FUNCTIONS ---------------------------------------------------- #
# Plots the upper bound set for a node and the lower bound set 
function plotBoundSets(UB::Vector{Constraint}, 
                       L::Vector{Solution{T}}) where T<:Real
    # Setup
    randNumber = rand(1:1000)
    println(randNumber)
    figure("Upper and lower bound sets"*string(randNumber),figsize=(6.5,5))
    xlabel(L"z^1(x)")
    ylabel(L"z^2(x)")
    PyPlot.title("Upper and lower bound sets")

    # Show the upper bound set 
    y_UBS1 = [y.point[1] for y in UB] 
    y_UBS2 = [y.point[2] for y in UB]
    scatter(y_UBS1, y_UBS2, color="green", marker="+", label = "UB")
    plot(y_UBS1, y_UBS2, color="green", linewidth=0.75, marker="+",
        markersize=1.0, linestyle=":")

    # Show the lower bound set 
    y_LBS1 = [y.z[1] for y in L] 
    y_LBS2 = [y.z[2] for y in L]
    scatter(y_LBS1, y_LBS2, color="red", marker="+", label = "L")
    plot(y_LBS1, y_LBS2, color="red", linewidth=0.75, marker="+",
        markersize=1.0, linestyle=":")

    #= display segments joining non-dominated points and their corners points
    Env1,Env2 = computeCornerPointsLowerEnvelop(y_LBS1, y_LBS2)
    plot(Env1, Env2, color="black", linewidth=0.75, marker="+", markersize=1.0, 
        linestyle=":")=#

    legend(bbox_to_anchor=[1,1], loc=0, borderaxespad=0, fontsize = "x-small")
end 

# Plot the nondominated points obtained by vOpt and the branch-and-bound method 
function plotYN(fname::String,
                ref::Vector{Vector{Float64}}, 
                L::Vector{Solution{T}}) where T<:Real
    # Setup
    figure("Nondominated points "*fname,figsize=(6.5,5))
    xlabel(L"z^1(x)")
    ylabel(L"z^2(x)")
    PyPlot.title("Nondominated points")

    # Show the set of nondominated points computed by vOpt 
    y_ref1 = [y[1] for y in ref] 
    y_ref2 = [y[2] for y in ref]
    scatter(y_ref1, y_ref2, color="black", marker="*", label = "ref")
    plot(y_ref1, y_ref2, color="black", linewidth=0.75, marker="*",
        markersize=1.0, linestyle=":")

    # Show the points computed by branch-and-bound
    y_N1 = [y.z[1] for y in L] 
    y_N2 = [y.z[2] for y in L]
    scatter(y_N1, y_N2, color="cyan", marker="+", label = "branch-and-bound")
    plot(y_N1, y_N2, color="cyan", linewidth=0.75, marker="+",
        markersize=1.0, linestyle=":")

    legend(bbox_to_anchor=[1,1], loc=0, borderaxespad=0, fontsize = "x-small")
end 

# Plot two upper bound sets
function plotUBS(ubdicho::Vector{Vector{Float64}}, 
                 ubparam::Vector{Vector{Float64}})

    fig = figure("Comparaison EBS dichotomique vs paramétrique", figsize=(6.5,5))
    plot([y[1] for y in ubdicho], [y[2] for y in ubdicho], color="green",
        linewidth=.75, marker="+", markersize=10., linestyle=":")
    plot([y[1] for y in ubparam], [y[2] for y in ubparam], color="red",
        linewidth=.75, marker="+", markersize=10., linestyle=":")
end 
        
# ----- VERIFICATION --------------------------------------------------------- #
using DataStructures

# Verifies that the upper bound set obtained by the parametric method is correct
# by comparing it to that obtained by the dichotomic method 
function verifyUBS(prob::_MOMKP, setvar::SetVariables, UB::Vector{Constraint})

    n::Int = size(prob.P)[2]
    nbSetVariables = length(setvar.setToOne) + length(setvar.setToZero)
    println("\nNumber of set variables : ", nbSetVariables)

    # Sum of the objective function values for the variables that have been
    # set to one 
    sumZ::Vector{Float64} = [0.,0.]

    # Residual capacity 
    ω_::Int = prob.ω[1]
    
    # Create a new variable containing the subproblem 
    newP = Matrix{Int}(undef, 2, n - nbSetVariables)
    newW = Matrix{Int}(undef, 1, n - nbSetVariables)

    inser = 1 
    for i in 1:n 
        if !(i in setvar.setToOne || i in setvar.setToZero)
            # The variable is added to the new subproblem
            newP[1,inser] = prob.P[1,i]
            newP[2,inser] = prob.P[2,i]
            newW[1,inser] = prob.W[1,i]
            inser += 1 

        elseif i in setvar.setToOne 
            sumZ += prob.P[:,i]
            ω_   -= prob.W[1,i]
        end 
    end 

    subProb = _MOMKP(newP, newW, [ω_])
    println(subProb)
    println("setvar : ", setvar)

    # Compute the upper bound set with the dichotomic method 
    init = initialisation(subProb, DICHOTOMIC)
    L    = Vector{Solution{Rational{Int}}}() 

    UBdicho::Vector{Constraint} = dichotomicMethod(subProb, init)

    paramPoints = OrderedSet([c.point for c in UB])
    dichoPoints = OrderedSet(sort([c.point + sumZ for c in UBdicho],rev=true))

    if paramPoints != dichoPoints 
        paramPoints = [y for y in paramPoints]
        dichoPoints = [y for y in dichoPoints]
        println("Parametric : ", paramPoints)
        println("Dichotomic : ", dichoPoints, "\n")
        plotUBS(dichoPoints, paramPoints)
    end 
end 

# Verifies that the initial solution and the lists of set variables are coherent
function verifySetvar(prob::_MOMKP,
                      setvar::SetVariables,
                      solInit::Solution{Float64})

    n = size(prob.P)[2]

    X  = zeros(Rational{Int},n)
    z  = [0.,0.]
    ω_ = prob.ω[1] 

    for i in 1:n 
        if i in setvar.setToOne
            X[i] = 1//1 
            z   += prob.P[:,i]
            ω_  -= prob.W[1,i]
        end 
    end 

    return X == solInit.X && z == solInit.z && ω_ == solInit.ω_ 
end 

# Prints whether the instance contains equivalent items for each instance in 
# the given directory 
function containsEquivalentItems(dir::String)

    files = readdir(dir*"dat/")

	for fname in files 

		println("\n", fname)

		# Read the instance in the file 
		if fname[length(fname)-3:length(fname)] == ".DAT"
			prob = readInstanceMOMKPformatPG(false, dir*"dat/"*fname)

		elseif fname[length(fname)-3:length(fname)] == ".dat"
			prob = readInstanceKP(dir*"dat/"*fname)

		else
			prob = readInstanceMOMKPformatZL(false, dir*"dat/"*fname)
		end

		newProb = groupEquivalentItems(prob)

        n1 = size(prob.P)[2]
        n2 = size(newProb.P)[2] 
        println("Contains equivalent items : ", !(n1 == n2))
	end 
end  

