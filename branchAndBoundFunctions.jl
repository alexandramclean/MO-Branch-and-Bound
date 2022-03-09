################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Branch-and-bound auxiliary functions                                 #
################################################################################

include("functions.jl")

@enum Order INCREASING DECREASING 

# ----- BRANCHING STRATEGIES ------------------------------------------------- #
# Returns the ranks of items ordered by decreasing order of utility  
function rank(u1::Vector{Rational{Int}}, # Utilities for objective function 1
              u2::Vector{Rational{Int}}) # Utilities for objective function 2 

    # Ranks for the first objective function 
    r1 = sortperm(sortperm(u1, rev=true))

    # Ranks for the second objective function 
    r2 = sortperm(sortperm(u2, rev=true))

    return r1, r2 
end 

# Items selected in increasing or decreasing order of min(u1,u2)
function minUtility(u1::Vector{Rational{Int}}, 
                    u2::Vector{Rational{Int}},
                    order::Order)

    min_utility = [min(u1[j], u2[j]) for j in 1:length(u1)]
    if order == DECREASING
        return sortperm(min_utility, rev=true) # Best variable first
    else 
        return sortperm(min_utility)           # Worst variable first
    end
end 

# Items selected in increasing or decreasing order of max(u1,u2) 
function maxUtility(u1::Vector{Rational{Int}}, 
                    u2::Vector{Rational{Int}},
                    order::Order)

    max_utility = [max(u1[j], u2[j]) for j in 1:length(u1)]
    if order == DECREASING 
        return sortperm(max_utility, rev=true) # Best variable first
    else 
        return sortperm(max_utility)           # Worst variable first

    end
end 

# Items selected in increasing or decreasing order of (u1+u2)/2 
function avgUtility(u1::Vector{Rational{Int}}, 
                    u2::Vector{Rational{Int}},
                    order::Order)

    avg_utility = [(u1[j] + u2[j])//2 for j in 1:length(u1)] 
    if order == DECREASING 
        return sortperm(avg_utility, rev=true) # Best variable first
    else 
        return sortperm(avg_utility)           # Worst variable first
    end 
end 

# Items selected in increasing or decreasing order of r1+r2 
function sumRank(r1::Vector{Int}, 
                 r2::Vector{Int},
                 order::Order) 
    
    sum_rank = [r1[j] + r2[j] for j in 1:length(r1)]
    if order == INCREASING 
        return sortperm(sum_rank)           # Best variable first
    else 
        return sortperm(sum_rank, rev=true) # Worst variable first
    end 
end 

# Items selected in increasing or decreasing order of min(r1,r2)
function minRank(r1::Vector{Int},
                 r2::Vector{Int},
                 order::Order)

    n = length(r1)
    min_rank = [min(r1[j],r2[j]) + (r1[j] + r2[j])//2*n for j in 1:n]
    if order == INCREASING 
        return sortperm(min_rank)           # Best variable first
    else
        return sortperm(min_rank, rev=true) # Worst variable first
    end 
end 

# Items selected in increasing or decreasing order of max(r1,r2)
function maxRank(r1::Vector{Int},
                 r2::Vector{Int},
                 order::Order)

    n = length(r1)
    max_rank = [max(r1[j],r2[j]) + (r1[j] + r2[j])//2*n for j in 1:n]
    if order == INCREASING 
        return sortperm(max_rank)           # Best variable first
    else
        return sortperm(max_rank, rev=true) # Worst variable first
    end 
end 

# ----- LOCAL NADIR POINTS AND DOMINANCE TESTS ------------------------------- # 
# Exemple
upperBound = DualBoundSet([[1.,8.], [4.,7.], [7.,4.], [9.,1.]], 
                           [Constraint(1//1, [9.,1.]), 
                           Constraint(7//15, [4.,7.]), 
                           Constraint(3//5, [7.,4.]),
                           Constraint(0//1, [1.,8.])])
incumbentSet1 = [[10.,4.], [9.,6.], [6.,8.], [3.,9.]] 
incumbentSet2 = [[10.,3.], [8.,5.], [4.,9.]]

# Computes the local nadir points 
function localNadirPoints(incumbentSet::Vector{Vector{Float64}})

    nadirs = Vector{Vector{Float64}}(undef, length(incumbentSet)-1)
    for i in 1:length(incumbentSet)-1 
        yl = incumbentSet[i+1]
        yr = incumbentSet[i] 
        nadirs[i] = [yl[1], yr[2]]
    end 
    return nadirs
end 

# Calcule les nadirs locaux décalés (hypothèse d'intégrité)
function shiftedLocalNadirPoints(nadirs::Vector)
    shiftedNadirs = Vector{Vector{Float64}}(undef, length(nadirs))
    for i in 1:length(nadirs)
        shiftedNadirs[i] = nadirs[i] + [1.,1.]
    end 
    return shiftedNadirs
end 

# Returns true if the local nadir point verifies all the constraints 
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

# Returns true if the upper bound set is dominated by the lower bound set 
function isDominated(UB::DualBoundSet, L::Vector{Vector{Float64}})

    is_dominated  = true 
    shiftedNadirs = shiftedLocalNadirPoints(localNadirPoints(L))
    i = 1 

    while is_dominated && i <= length(shiftedNadirs)
        is_dominated = is_dominated && verifiesConstraints(UB.constraints, shiftedNadirs[i])
        i += 1 
    end 
    return is_dominated
end 
