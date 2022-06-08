################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Data structures                                                      #
################################################################################

@enum Optimisation MAX MIN

@enum Method PARAMETRIC_LP DICHOTOMIC SIMPLEX PARAMETRIC_MT 

# ----- PROBLEMS ------------------------------------------------------------- #
# Datastructure of a multi-objective multi-dimensionnal KP with 0/1 variables
struct _MOMKP
    P  :: Matrix{Int} # profit of items for the objectives, k=1..p, j=1..n
    W  :: Matrix{Int} # weight of items for the constraints, i=1..m, j=1..n
    ω  :: Vector{Int} # capacity of knapsacks, i=1..m
end

function didacticInstance()
    return _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11]) 
end 

function example()
    return _MOMKP([11 2 2 8 10 9 1 16 ; 2 7 7 8 4 1 3 4], [4 4 4 6 4 3 2 6], [13])
end 

# ----- SOLUTIONS ------------------------------------------------------------ #
# Data structure of a solution that can contain fractions of an object
mutable struct Solution{T<:Real} 
    X::Vector{Rational{Int}} # Vector of binary variables 
    z::Vector{T}             # Values for the objective functions 
    ω_::Int                  # Residual capacity 
end 

# Creates an initial solution for prob with no set variables 
Solution{Float64}(prob::_MOMKP) = 
    Solution(zeros(Rational{Int}, size(prob.P)[2]), [0.,0.], prob.ω[1])
Solution{Rational{Int}}(prob::_MOMKP) = 
    Solution(zeros(Rational{Int}, size(prob.P)[2]), [0//1,0//1], prob.ω[1])
                                
Solution(t::Vector{Float64}) = Solution([0//1], t, 0)

# Transform a solution of rational type to a solution of float type 
Solution{Float64}(x::Solution{Rational{Int}}) = 
    Solution(x.X, [Float64(i) for i in x.z], x.ω_)

# ----- TRANSPOSITIONS ------------------------------------------------------- #
# Stores a critical weight λ and the corresponding transposition(s)
struct Transposition
	λ::Rational{Int}			  # Critical weight 
	pairs::Vector{Tuple{Int,Int}} # List of pairs of items to swap 
end

Transposition(λ) = Transposition(λ,[])

# Stores the transpositions and initial sequence and positions 
struct Initialisation
    r1::Union{Vector{Rational{Int}},Nothing} # Utilities for z1 
    r2::Union{Vector{Rational{Int}},Nothing} # Utilities for z2
    # Critical weights and transpositions
    transpositions::Union{Vector{Transposition},Nothing} 
    seq::Union{Vector{Int},Nothing}          # Initial sequence 
    pos::Union{Vector{Int},Nothing}          # Positions in the sequence
end 

# Structure containing the list of variables that have been set as well as 
# the index of the remaining transpositions 
struct SetVariables 
    setToOne::Vector{Int}          # List of variables set to one 
    setToZero::Vector{Int}         # List of variables set to zero 
    # Index of remaining transpositions + sequence and positions 
    # (for the parametric method)
    transpInd::Union{Vector{Vector{Int}},Nothing}  
    seq::Union{Vector{Int},Nothing} # Initial sequence with set variables
    pos::Union{Vector{Int},Nothing} # Positions of the objects in the sequence 
end 

# ----- BOUND SETS ----------------------------------------------------------- #
# Data structure of a constraint generated while computing the upper bound set
struct Constraint 
    # Critical weight
	λ::Union{Rational{Int},Float64} 
    # Associated point (u1,u2)
	point::Union{Vector{Rational{Int}},Vector{Float64}} 
end
# The associated constraint is λz1 + (1-λ)z2 <= λu1 + (1-λ)u2

# Data structure representing the dual bound set (upper bound set in this case)  
struct DualBoundSet{T} 
    points::Vector{Vector{T}}    
    constraints::Vector{Constraint}
end

DualBoundSet{Float64}()       = DualBoundSet(Vector{Float64}[], Constraint[]) 
DualBoundSet{Rational{Int}}() = DualBoundSet(Vector{Rational{Int}}[], Constraint[])


# ----- BRANCH-AND-BOUND ----------------------------------------------------- #
@enum Status DOMINANCE OPTIMALITY INFEASIBILITY NOTPRUNED MAXDEPTH

# Data structure representing a node in a branch-and-bound algorithm
mutable struct Node 
    UB::Union{DualBoundSet,Nothing} # Upper bound set for the node
    setvar::Union{SetVariables,Initialisation} # Initialisation for the node 
    solInit::Solution               # Initial solution with set variables
    # Indicates whether the node has been pruned and for what reason 
    status::Status                   
end