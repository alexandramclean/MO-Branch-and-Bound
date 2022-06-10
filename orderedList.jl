################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Liste ordonnée de points non-dominés                                 #
################################################################################

using Random

# ----- AUXILIARY FUNCTIONS -------------------------------------------------- #
# Precondition : x and y are of the same type 
# Returns true if x is strictly smaller than y on dimension dim 
function isStrictlySmaller(x::Union{Solution{T}, Vector{T}},
                           y::Union{Solution{T}, Vector{T}},
                           dim::Int) where T<:Real
    if typeof(x) == Solution{Float64} || typeof(x) == Solution{Rational{Int}} 
        return x.z[dim] < y.z[dim]
    else 
        return x[dim] < y[dim] 
    end 
end 

# ----- VERIFY --------------------------------------------------------------- #
# Searches for the index of the last point dominated by y 
function lastDominatedPoint(yN::Vector{Vector{T}},
                            y::Vector{T},
                            start::Int64, 
                            finish::Int64,
                            opt::Optimisation=MAX) where T<:Real
    if start >= finish
        return start
    elseif opt == MIN && dominates(y, yN[finish], opt)
        return finish
    elseif opt == MAX && dominates(y, yN[start], opt)
        return start 
    end
    
    mid = div(start + finish, 2)
    
    if dominates(y, yN[mid], opt)
        if opt == MIN
            if !dominates(y, yN[mid+1], opt)
                # yN[mid] is the last dominated point 
                return mid
            else 
                return lastDominatedPoint(yN, y, mid+1, finish, opt)
            end 
        elseif opt == MAX
            if !dominates(y, yN[mid-1], opt)
                # yN[mid] is the first dominated point 
                return mid 
            else 
                return lastDominatedPoint(yN, y, start, mid-1, opt) 
            end 
        end
        
    # y does not dominate yN[mid]
    elseif opt == MIN
        return lastDominatedPoint(yN, y, start, mid-1, opt)
    elseif opt == MAX 
        return lastDominatedPoint(yN, y, mid+1, finish, opt)
    end
end

# Verifies that y's successors (>ind if MIN or <ind if MAX) are not dominated 
function verify(yN::Union{Vector{Solution{T}}, Vector{Vector{T}}},
                y::Union{Solution{T}, Vector{T}}, 
                ind::Int64,
                opt::Optimisation=MAX) where T<:Real
                
    indLastDominated = 0

    # If y dominates any points it dominates its immediate successor 
    # (or predecessor if opt == MAX)

    if opt == MIN && ind < length(yN) && dominates(y, yN[ind+1], opt) 

        # The ordered list contains solutions 
        if typeof(y) == Solution{Float64} || typeof(y) == Solution{Rational{Int}}
            indLastDominated = lastDominatedPoint([sol.z for sol in yN], y.z, ind+1, length(yN), opt)
        else 
            # The ordered list contains points 
            indLastDominated = lastDominatedPoint(yN, y, ind+1, length(yN), opt)
        end 

        # Deletion of the dominated points 
        for j in indLastDominated:-1:ind+1
            deleteat!(yN, j)
        end

    elseif opt == MAX && ind > 1 && dominates(y, yN[ind-1], opt)
        # First dominated point

        # The ordered list contains solutions 
        if typeof(y) == Solution{Float64} || typeof(y) == Solution{Rational{Int}}
            indLastDominated = lastDominatedPoint([sol.z for sol in yN], y.z, 1, ind-1, opt)
        else 
            # The ordered list contains points 
            indLastDominated = lastDominatedPoint(yN, y, 1, ind-1, opt)
        end 

        # Deletion of the dominated points 
        for j in ind-1:-1:indLastDominated
            deleteat!(yN, j)
        end 
    end 
end

# ----- ADD ------------------------------------------------------------------ #
# Returns the index of y so as to maintain an ordered list 
# If the index returned is 0 then y should not be inserted (dominated point or 
# duplicate)
function addRec(yN::Vector{Vector{T}},
                y::Vector{T},
                start::Int64, 
                finish::Int64,
                opt::Optimisation=MAX) where T<:Real
                     
    # Stopping criterion : only 1 element left in the list            
    if start >= finish 

        # y is inserted before yN[start] and yN[start] will be deleted
        if dominates(y, yN[start], opt)
            if opt == MIN
                return start
            else 
                return start+1
            end 
        
        # y is dominated and not inserted
        elseif dominates(yN[start], y, opt)
            return 0 

        # There is no dominance between y and yN[start]
        elseif y[1] < yN[start][1] 
            #!isStrictlySmaller(y, yN[start], 1)
            # y is inserted before yN[start]
            return start 

        else # y is inserted after yN[start]
            return start+1
        end
    end
    
    mid = div(start+finish,2)
    
    # If y is dominated by a point in the list it is not inserted
    if (dominates(yN[mid], y, opt) 
        || dominates(yN[start], y, opt) 
        || dominates(yN[finish], y, opt))
        return 0

    elseif y[1] < yN[mid][1] 
        #!isStrictlySmaller(y, yN[mid], 1)
        # The case where y and yN[mid] have equal values for both objective 
        # functions is included in the dominance test 
        return addRec(yN, y, start, mid-1, opt)
    else
        return addRec(yN, y, mid+1, finish, opt)
    end
end

# Adds an element to an ordered list 
# The ordered list can be an array of arrays (points) or an array 
# of solutions
function add!(yN::Union{Vector{Solution{T}}, Vector{Vector{T}}}, 
              y::Union{Solution{T},Vector{T}}, 
              opt::Optimisation=MAX) where T<:Real

    # The ordered list contains solutions 
    if typeof(y) == Solution{Float64} || typeof(y) == Solution{Rational{Int}}
        #! Should not be needed 
        y.z = [floor(y.z[1]), floor(y.z[2])]

        # Search for the position of y 
        if length(yN) == 0
            insert!(yN, 1, y)

        else
            # addRecc only accepts an array of points 
            ind = addRec([sol.z for sol in yN], y.z, 1, length(yN), opt)
            if ind > 0 
                insert!(yN, ind, y)

                # Elimination of the dominated points 
                verify(yN, y, ind, opt)
            end 
        end 
    else 
        # The ordered list contains points 
        #! Should not be required 
        if typeof(y) == Vector{Float64}
            y = [floor(y[1]), floor(y[2])]
        end 

        # Search for the position of y 
        if length(yN) == 0
            insert!(yN, 1, y)

        else
            # addRecc only accepts an array of points 
            ind = addRec(yN, y, 1, length(yN), opt)
            if ind > 0 
                insert!(yN, ind, y)

                # Elimination of the dominated points 
                verify(yN, y, ind, opt)
            end 
        end 
    end 
end

# ----- PRINT ---------------------------------------------------------------- #
# Prints yN in a readable format
function afficher(yN::Vector{Vector{T}}) where T<:Real
    print("| ")
    for i in 1:length(yN)
        if yN[i][1] < 10
            print("  ", yN[i][1], "  ")
        elseif yN[i][1] < 100
            print(" ", yN[i][1], "  ")
        else
            print(yN[i][1], "  ")
        end
    end
    print("|\n| ")
    for i in 1:length(yN)
        if yN[i][2] < 10
            print("  ", yN[i][2], "  ")
        elseif yN[i][2] < 100
            print(" ", yN[i][2], "  ")
        else
            print(yN[i][2], "  ")
        end
    end
    print("|\n")
end
