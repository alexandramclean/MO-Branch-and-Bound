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
# Computes the local nadir points 
function localNadirPoints(incumbentSet::Vector{Solution})

    nadirs = Vector{Vector{Float64}}(undef, length(incumbentSet)-1)
    for i in 1:length(incumbentSet)-1 
        yl = incumbentSet[i].z
        yr = incumbentSet[i+1].z 
        nadirs[i] = [yl[1], yr[2]]
    end 
    return nadirs
end 

# Calcule les nadirs locaux décalés (hypothèse d'intégrité)
function shiftedLocalNadirPoints(nadirs::Vector{Vector{Float64}})

    shiftedNadirs = Vector{Vector{Float64}}(undef, length(nadirs))
    for i in 1:length(nadirs)
        shiftedNadirs[i] = nadirs[i] + [1.,1.]
    end 
    return shiftedNadirs
end 

# Returns true if the point y verifies all the constraints 
function verifiesConstraints(constraints::Vector{Constraint}, 
                             y::Vector{Float64})
            
    verif = true 
    i = 1 

    while verif && i <= length(constraints)
        λ     = constraints[i].λ 
        point = constraints[i].point
        verif = verif && λ*y[1] + (1-λ)*y[2] <= λ*point[1] + (1-λ)*point[2]
        i    += 1 
    end 
    return verif 
end 

# Plots the upper bound set for a node and the lower bound set 
function plotBoundSets(UB::DualBoundSet, L::Vector{Solution})
    # Setup
    randNumber = rand(1:100)
    println(randNumber)
    figure("Upper and lower bound sets"*string(randNumber),figsize=(6.5,5))
    xlabel(L"z^1(x)")
    ylabel(L"z^2(x)")
    PyPlot.title("Upper and lower bound sets")

    # Show the upper bound set 
    y_UBS1 = [y[1] for y in UB.points] 
    y_UBS2 = [y[2] for y in UB.points]
    scatter(y_UBS1, y_UBS2, color="green", marker="+", label = "UB")
    plot(y_UBS1, y_UBS2, color="green", linewidth=0.75, marker="+",
        markersize=1.0, linestyle=":")

    # Show the points computed by the dichotomic method 
    y_LBS1 = [y.z[1] for y in L] 
    y_LBS2 = [y.z[2] for y in L]
    scatter(y_LBS1, y_LBS2, color="red", marker="+", label = "L")
    plot(y_LBS1, y_LBS2, color="red", linewidth=0.75, marker="+",
        markersize=1.0, linestyle=":")
end 


# Returns true if the upper bound set for a particular node is dominated by the 
# lower bound set 
function isDominated(UB::DualBoundSet, L::Vector{Solution})

    is_dominated  = true 
    shiftedNadirs = shiftedLocalNadirPoints(localNadirPoints(L))
    i = 1 

    while is_dominated && i <= length(shiftedNadirs)
        is_dominated = is_dominated && !verifiesConstraints(UB.constraints, shiftedNadirs[i])
        i += 1 
    end 
    return is_dominated
end 
