################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Liste ordonnée de points non-dominés                                 #
################################################################################

using Random

# ----- AUXILIARY FUNCTIONS -------------------------------------------------- #
# Precondition : x and y are of the same type 
# Returns true if x is strictly smaller than y on dimension dim 
function isStrictlySmaller(x::Union{Solution, Vector{T}},
                           y::Union{Solution, Vector{T}},
                           dim::Int) where T<:Real
    if typeof(x) == Solution{Float64} || typeof(x) == Solution{Rational{Int}} 
        return x.z[dim] < y.z[dim]
    else 
        return x[dim] < y[dim] 
    end 
end 

# ----- VERIFY --------------------------------------------------------------- #
# Searches for the index of the last point dominated by y 
function lastDominatedPoint(yN::Union{Vector{Solution}, Vector{Vector{T}}},
                            y::Union{Solution, Vector{T}},
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
        if ((opt == MIN && !dominates(y, yN[mid+1], opt))
            || (opt == MAX && !dominates(y, yN[mid-1], opt)))
            # yN[mid] is the last dominated point 
            return mid
        elseif opt == MIN
            return lastDominatedPoint(yN, y, mid+1, finish, opt)
        else
            return lastDominatedPoint(yN, y, start, mid-1, opt)
        end
        
    # y does not dominate yN[mid]
    elseif opt == MIN
        return lastDominatedPoint(yN, y, start, mid-1, opt)
    else
        return lastDominatedPoint(yN, y, mid+1, finish, opt)
    end
end

# Verifies that y's successors (>ind) are not dominated by y 
# Vérification que les successeurs de y (>ind) ne sont pas dominés par y
# Si y domine des points de yN alors il domine sont successeur large minimum
function verify(yN::Union{Vector{Solution}, Vector{Vector{T}}},
                y::Union{Solution, Vector{T}}, 
                ind::Int64,
                opt::Optimisation=MAX) where T<:Real
                  
    if opt == MIN && ind < length(yN) && dominates(y, yN[ind+1], opt)
        indLastDominated = lastDominatedPoint(yN, y, ind+1, length(yN), opt)
        #println("Last dominated point : ", indLastDominated)
        for j in indLastDominated:-1:ind+1
            deleteat!(yN, j)
        end
        
    elseif opt == MAX && ind > 1 && dominates(y, yN[ind-1], opt)
        indFirstDominated = lastDominatedPoint(yN, y, 1, ind-1, opt)
        #println("First dominated point : ", indFirstDominated)
        for j in ind-1:-1:indFirstDominated
            deleteat!(yN, j)
        end
    end
end

# ----- ADD ------------------------------------------------------------------ #
function addRec(yN::Union{Vector{Solution}, Vector{Vector{T}}},
                y::Union{Solution, Vector{T}},
                start::Int64, 
                finish::Int64,
                opt::Optimisation=MAX) where T<:Real
                     
    # Stopping criterion : only 1 element left in the list            
    if start >= finish 
        if dominates(yN[start], y, opt)
            return 0
        elseif isStrictlySmaller(yN[start], y, 1)
            return start + 1
        else
            return start
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
        # Le cas d'égalité sur les deux fonctions objectif a déjà été traité
        # dans le test de dominance
        return addRec(yN, y, start, mid-1, opt)
    else
        return addRec(yN, y, mid+1, finish, opt)
    end
end

function add!(yN::Union{Vector{Solution}, Vector{Vector{T}}}, 
              y::Union{Solution,Vector{T}}, 
              opt::Optimisation=MAX) where T<:Real
    #! Affichage
    #println("Adding ", y)

    # Search for the position of y 
    if length(yN) == 0
        insert!(yN, 1, y)
        #afficher(yN)
    else
        ind = addRec(yN, y, 1, length(yN), opt)
        if ind > 0 
            insert!(yN, ind, y)
            #afficher(yN)

            # Elimination of the dominated points
            verify(yN, y, ind, opt)
            #println("Verification : ")
            #afficher(yN)
        end
    end
end

# ----- AFFICHER ------------------------------------------------------------- #
function afficher(yN)
    print("|  ")
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
end
