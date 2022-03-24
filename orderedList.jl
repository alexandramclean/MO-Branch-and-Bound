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
function lastDominatedPoint(yN::Union{Vector{Solution{T}}, Vector{Vector{T}}},
                            y::Union{Solution{T}, Vector{T}},
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

# Verifies that y's successors (>ind) are not dominated by y 
# Vérification que les successeurs de y (>ind) ne sont pas dominés par y
# Si y domine des points de yN alors il domine sont successeur large minimum
function verify(yN::Union{Vector{Solution{T}}, Vector{Vector{T}}},
                y::Union{Solution{T}, Vector{T}}, 
                ind::Int64,
                opt::Optimisation=MAX) where T<:Real
                
    indLastDominated = 0

    # If y dominates any points it dominates its immediate successor 
    # (or predecessor if opt == MAX)
    if opt == MIN && ind < length(yN) && dominates(y, yN[ind+1], opt) 
        indLastDominated = lastDominatedPoint(yN, y, ind+1, length(yN), opt)
        #println("Last dominated point : ", indLastDominated)
        for j in indLastDominated:-1:ind+1
            deleteat!(yN, j)
        end
    elseif opt == MAX && ind > 1 && dominates(y, yN[ind-1], opt)
        # First dominated point
        indLastDominated = lastDominatedPoint(yN, y, 1, ind-1, opt)

        for j in ind-1:-1:indLastDominated
            deleteat!(yN, j)
        end 
    end 

    return indLastDominated
end

# ----- LOCAL NADIR POINTS --------------------------------------------------- #
# Retourne un tableau contenant les indices des points nadirs locaux 
# correspondant au point yN.solutions[ind]
function correspondingNadirPoints(yN::PrimalBoundSet{T},
                                  ind::Int) where T<:Real 

    if ind == 1
        return [ind] 
    elseif ind == length(yN.solutions)
        return [ind-1]
    else 
        return [ind-1, ind]
    end 
end 

function nadirs!(yN::PrimalBoundSet{T}, # Lower bound set 
                 y::Solution{T},        # Solution that was added 
                 ind::Int,              # Index of the solution that was added 
                 indLastDominated::Int, # Index of the last (or first) point 
                 # dominated by y
                 opt::Optimisation) where T<:Real

    if indLastDominated == 0 # No solutions were deleted 

        if opt == MAX && ind < length(yN.solutions) && ind > 1 
            # The local nadir corresponding to the points before and after ind
            # needs to be deleted 
            deleteat!(yN.nadirs, ind-1)

            # Two new local nadir points 
            yl = yN.solutions[ind-1]
            yr = yN.solutions[ind+1]

            insert!(yN.nadirs, ind-1, [yl.z[1], y.z[2]])
            insert!(yN.nadirs, ind, [y.z[1], yr.z[2]])    

        elseif ind == 1
            # No local nadirs need to be deleted 
            # A new local nadir is inserted at the start of the list 
            yr = yN.solutions[ind+1]
            insert!(yN.nadirs, 1, [y.z[1], yr.z[2]])

        elseif ind == length(yN.solutions)  
            # No local nadirs need to be deleted 
            # A new local nadir is inserted at the end of the list 
            yl = yN.solutions[ind-1] 
            insert!(yN.nadirs, length(yN.nadirs)+1, [yl.z[1], y.z[2]])      
        end 
    else # Some solutions have been deleted 

        # The solutions between indLastDominated and ind-1 have been deleted
        # The corresponding local nadir points are deleted 
        if ind == length(yN.solutions)+1
            correspondingNadirs = indLastDominated-1:length(yN.solutions)-1
        elseif indLastDominated == 1
            correspondingNadirs = indLastDominated:ind-1 
        else 
            correspondingNadirs = indLastDominated-1:ind-1
        end 
        println(correspondingNadirs)

        for i in reverse(correspondingNadirs)
            deleteat!(yN.nadirs, i)
        end 

        # Add the new local nadir points 
        if indLastDominated == 1
            # A new local nadir is inserted at the start of the list 
            yr = yN.solutions[indLastDominated+1]
            insert!(yN.nadirs, 1, [y.z[1], yr.z[2]])

        elseif indLastDominated == length(yN.solutions)
            # A new local nadir is inserted at the end of the list
            yl = yN.solutions[indLastDominated-1]
            insert!(yN.nadirs, length(yN.nadirs)+1, [yl.z[1], y.z[2]])
            
        else 
            # Two new local nadir points are inserted 
            yl = yN.solutions[indLastDominated-1]
            yr = yN.solutions[indLastDominated+1]

            insert!(yN.nadirs, indLastDominated, [yl.z[1], y.z[2]])
            insert!(yN.nadirs, indLastDominated+1, [y.z[1], yr.z[2]])
        end 
    end 
end 

# ----- ADD ------------------------------------------------------------------ #
function addRec(yN::Union{Vector{Solution{T}}, Vector{Vector{T}}},
                y::Union{Solution{T}, Vector{T}},
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
        elseif isStrictlySmaller(y, yN[start], 1)
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

    elseif isStrictlySmaller(y, yN[mid], 1)
        # The case where y and yN[mid] have equal values for both objective 
        # functions is included in the dominance test 
        return addRec(yN, y, start, mid-1, opt)
    else
        return addRec(yN, y, mid+1, finish, opt)
    end
end

function add!(yN::Union{PrimalBoundSet{T}, Vector{Vector{T}}}, 
              y::Union{Solution{T},Vector{T}}, 
              opt::Optimisation=MAX) where T<:Real

    if typeof(yN) == PrimalBoundSet{Float64} || 
        typeof(yN) == PrimalBoundSet{Rational{Int}}

        #! Affichage
        println("\nAdding ", y.z)

        # Search for the position of y 
        if length(yN.solutions) == 0
            insert!(yN.solutions, 1, y)

        else
            ind = addRec(yN.solutions, y, 1, length(yN.solutions), opt)
            println("ind = ", ind)
            if ind > 0 
                insert!(yN.solutions, ind, y)

                # Elimination of the dominated points 
                indLastDominated = verify(yN.solutions, y, ind, opt)
                println("indLastDominated = ", indLastDominated)

                # Recalculer les nadirs locaux affectés
                if length(yN.nadirs) == 0
                    # Initialiser la liste des nadirs locaux 
                    yN.nadirs = localNadirPoints(yN.solutions)
                elseif opt == MIN 
                    nadirs!(yN, y, ind, indLastDominated, opt)
                else 
                    nadirs!(yN, y, ind, indLastDominated, opt)
                end 
            end 
        end 
    else 
        # Search for the position of y 
        if length(yN) == 0
            insert!(yN, 1, y)
            #afficher(yN)
        else
            #println("length = ", length(yN))
            ind = addRec(yN, y, 1, length(yN), opt)
            #println("ind = ", ind)
            if ind > 0 
                insert!(yN, ind, y)
                #afficher(yN)

                # Elimination of the dominated points
                #println("Verification : ")
                verify(yN, y, ind, opt)
            end
            #afficher(yN)
        end
    end
end

# ----- AFFICHER ------------------------------------------------------------- #
# Prints yN in a readable format
function afficher(yN)
    if typeof(yN[1]) == Solution{Float64} || 
            typeof(yN[1]) == Solution{Rational{Int}}
        print("| ")
        for i in 1:length(yN)
            if yN[i].z[1] < 10
                print("  ", yN[i].z[1], "  ")
            elseif yN[i].z[1] < 100
                print(" ", yN[i].z[1], "  ")
            else
                print(yN[i].z[1], "  ")
            end
        end
        print("|\n| ")
        for i in 1:length(yN)
            if yN[i].z[2] < 10
                print("  ", yN[i].z[2], "  ")
            elseif yN[i].z[2] < 100
                print(" ", yN[i].z[2], "  ")
            else
                print(yN[i].z[2], "  ")
            end
        end
        print("|\n")
    else 
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
end

# Exemple 
L = PrimalBoundSet{Float64}(
    [Solution{Float64}(Rational{Int}[], [3., 9.], 0),
     Solution{Float64}(Rational{Int}[], [6., 8.], 0),
     Solution{Float64}(Rational{Int}[], [9., 6.], 0),
     Solution{Float64}(Rational{Int}[], [10., 4.], 0)
    ], [[3.,8.], [6.,6.], [9.,4.]])
