################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Liste ordonnée de points non-dominés                                 #
################################################################################

include("dataStructures.jl")
using Random
@enum Optimisation Max Min

# Retourne vrai si x domine y
function domine(x, y, opt::Optimisation=Max)
    if opt == Min
        return ((x.z[1] <= y.z[1] && x.z[2] < y.z[2])
            || (x.z[1] < y.z[1] && x.z[2] <= y.z[2])
            || (x.z[1] == y.z[1] && x.z[2] == y.z[2])) # Pas de doublons
    else
        return ((x.z[1] >= y.z[1] && x.z[2] > y.z[2])
            || (x.z[1] > y.z[1] && x.z[2] >= y.z[2])
            || (x.z[1] == y.z[1] && x.z[2] == y.z[2])) # Pas de doublons
    end
end

## VERIFIER
# Recherche de l'indice du dernier point dominé par y
function dernier_point_domine(sN::Vector{Solution},
                              y::Solution,
                              deb::Int64, fin::Int64,
                              opt::Optimisation=Max) 
    if deb >= fin
        return deb
    elseif opt == Min && domine(y, sN[fin], opt)
        return fin
    elseif opt == Max && domine(y, sN[deb], opt)
        return deb
    end
    
    mil = div(deb + fin, 2)
    
    if domine(y, sN[mil], opt)
        if ((opt == Min && !domine(y, sN[mil+1], opt))
            || (opt == Max && !domine(y, sN[mil-1], opt)))
            # sN[mil] est le dernier point dominé
            return mil
        elseif opt == Min
            return dernier_point_domine(sN, y, mil+1, fin, opt)
        else
            return dernier_point_domine(sN, y, deb, mil-1, opt)
        end
        
    # y ne domine pas sN[mil]
    elseif opt == Min
        return dernier_point_domine(sN, y, deb, mil-1, opt)
    else
        return dernier_point_domine(sN, y, mil+1, fin, opt)
    end
end

# Vérification que les successeurs de y (>ind) ne sont pas dominés par y
# Si y domine des points de sN alors il domine sont successeur large minimum
function verifier(sN::Vector{Solution},
                  y::Solution, ind::Int64,
                  opt::Optimisation=Max) 
                  
    if opt == Min && ind < length(sN) && domine(y, sN[ind+1], opt)
        ind_dernier_domine = dernier_point_domine(sN, y, ind+1, length(sN), opt)
        #println("Dernier dominé : ", ind_dernier_domine)
        for j in ind_dernier_domine:-1:ind+1
            deleteat!(sN, j)
        end
        
    elseif opt == Max && ind > 1 && domine(y, sN[ind-1], opt)
        ind_premier_domine = dernier_point_domine(sN, y, 1, ind-1, opt)
        #println("Premier dominé : ", ind_premier_domine)
        for j in ind-1:-1:ind_premier_domine
            deleteat!(sN, j)
        end
    end
end

## AJOUTER
function ajouter_rec(sN::Vector{Solution},
                     y::Solution,
                     deb::Int64, fin::Int64,
                     opt::Optimisation=Max)
                     
    if deb >= fin # Cas d'arrêt : 1 seul élément dans la sous-liste
        if domine(sN[deb], y, opt)
            return 0
        elseif sN[deb].z[1] < y.z[1]
            return deb + 1
        else
            return deb
        end
    end
    
    mil = div(deb+fin,2)
    
    # Si y est dominé par un point dans la liste on ne l'insère pas
    if domine(sN[mil], y, opt) || domine(sN[deb], y, opt) || domine(sN[fin], y, opt)
        return 0
    elseif sN[mil].z[1] > y.z[1] #|| sN[mil].z[1] == y.z[1]
        # Le cas d'égalité sur les deux fonctions objectif a déjà été traité
        # dans le test de dominance
        return ajouter_rec(sN, y, deb, mil-1, opt)
    else
        return ajouter_rec(sN, y, mil+1, fin, opt)
    end
end

function ajouter(sN::Vector{Solution}, y::Solution, opt::Optimisation=Max)
    #! Affichage
    #println("Ajout de ", y.z)

    # Recherche de la position de y
    if length(sN) == 0
        insert!(sN, 1, y)
        #afficher(sN)
    else
        ind = ajouter_rec(sN, y, 1, length(sN), opt)
        if ind > 0 # Insertion
            insert!(sN, ind, y)
            #afficher(sN)

            # Elimination des points dominés
            verifier(sN, y, ind, opt)
            #println("Vérification : ")
            #afficher(sN)
        end
    end
end

## AFFICHER
function afficher(sN)
    print("|  ")
    for i in 1:length(sN)
        if sN[i].z[1] < 10
            print("  ", sN[i].z[1], "  ")
        elseif sN[i].z[1] < 100
            print(" ", sN[i].z[1], "  ")
        else
            print(sN[i].z[1], "  ")
        end
    end
    print("|\n| ")
    for i in 1:length(sN)
        if sN[i].z[2] < 10
            print("  ", sN[i].z[2], "  ")
        elseif sN[i].z[2] < 100
            print(" ", sN[i].z[2], "  ")
        else
            print(sN[i].z[2], "  ")
        end
    end
    print("|\n")
end
