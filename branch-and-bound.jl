################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Branch-and-bound algorithm                                           #
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
        return sortperm(min_utility, rev=true)
    else 
        return sortperm(min_utility)
    end
end 

# Items selected in increasing or decreasing order of max(u1,u2) 
function maxUtility(u1::Vector{Rational{Int}}, 
                    u2::Vector{Rational{Int}},
                    order::Order)

    max_utility = [max(u1[j], u2[j]) for j in 1:length(u1)]
    if order == DECREASING 
        return sortperm(max_utility, rev=true)
    else 
        return sortperm(max_utility)
end

# Items selected in increasing or decreasing order of (u1+u2)/2 
function avgUtility(u1::Vector{Rational{Int}}, 
                    u2::Vector{Rational{Int}},
                    order::Order)

    avg_utility = [(u1[j] + u2[j])//2 for j in 1:length(u1)] 
    if order == DECREASING 
        return sortperm(avg_utility, rev=true)
    else 
        return sortperm(avg_utility)
    end 
end 

# Items selected in increasing or decreasing order of r1+r2 
function sumRank(r1::Vector{Int}, 
                 r2::Vector{Int},
                 order::Order) 
    
    sum_rank = [r1[j] + r2[j] for j in 1:length(r1)]
    if order == DECREASING 
        return sortperm(sum_rank, rev=true) 
    else 
        return sortperm(sum_rank)
    end 
end 

# Items selected in increasing or decreasing order of min(r1,r2)
function minRank(r1::Vector{Int},
                 r2::Vector{Int},
                 order::Order)

    min_rank = [min(r1[j],r2[j]) for j in 1:length(r1)]
    if order == DECREASING 
        return sortperm(min_rank, rev=true)
    else
        return sortperm(min_rank)
    end 
end 

# Items selected in increasing or decreasing order of max(r1,r2)
