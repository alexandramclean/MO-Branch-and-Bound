################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Liste ordonnée de points non-dominés                                 #
################################################################################

using Random

## VERIFIER
# Recherche de l'indice du dernier point dominé par y
function dernier_point_domine(yN::Vector{Solution},
                              y::Solution,
                              deb::Int64, fin::Int64,
                              opt::Optimisation=MAX) 
    if deb >= fin
        return deb
    elseif opt == MIN && domine(y, yN[fin], opt)
        return fin
    elseif opt == MAX && domine(y, yN[deb], opt)
        return deb
    end
    
    mil = div(deb + fin, 2)
    
    if domine(y, yN[mil], opt)
        if ((opt == MIN && !domine(y, yN[mil+1], opt))
            || (opt == MAX && !domine(y, yN[mil-1], opt)))
            # yN[mil] est le dernier point dominé
            return mil
        elseif opt == MIN
            return dernier_point_domine(yN, y, mil+1, fin, opt)
        else
            return dernier_point_domine(yN, y, deb, mil-1, opt)
        end
        
    # y ne domine pas yN[mil]
    elseif opt == MIN
        return dernier_point_domine(yN, y, deb, mil-1, opt)
    else
        return dernier_point_domine(yN, y, mil+1, fin, opt)
    end
end

# Vérification que les successeurs de y (>ind) ne sont pas dominés par y
# Si y domine des points de yN alors il domine sont successeur large minimum
function verifier(yN::Vector{Solution},
                  y::Solution, 
                  ind::Int64,
                  opt::Optimisation=MAX) 
                  
    if opt == MIN && ind < length(yN) && domine(y, yN[ind+1], opt)
        ind_dernier_domine = dernier_point_domine(yN, y, ind+1, length(yN), opt)
        #println("Dernier dominé : ", ind_dernier_domine)
        for j in ind_dernier_domine:-1:ind+1
            deleteat!(yN, j)
        end
        
    elseif opt == MAX && ind > 1 && domine(y, yN[ind-1], opt)
        ind_premier_domine = dernier_point_domine(yN, y, 1, ind-1, opt)
        #println("Premier dominé : ", ind_premier_domine)
        for j in ind-1:-1:ind_premier_domine
            deleteat!(yN, j)
        end
    end
end

## AJOUTER
function ajouter_rec(yN::Vector{Solution},
                     y::Solution,
                     deb::Int64, fin::Int64,
                     opt::Optimisation=MAX)
                     
    if deb >= fin # Cas d'arrêt : 1 seul élément dans la sous-liste
        if domine(yN[deb], y, opt)
            return 0
        elseif yN[deb].z[1] < y.z[1]
            return deb + 1
        else
            return deb
        end
    end
    
    mil = div(deb+fin,2)
    
    # Si y est dominé par un point dans la liste on ne l'insère pas
    if domine(yN[mil], y, opt) || domine(yN[deb], y, opt) || domine(yN[fin], y, opt)
        return 0
    elseif yN[mil].z[1] > y.z[1] 
        # Le cas d'égalité sur les deux fonctions objectif a déjà été traité
        # dans le test de dominance
        return ajouter_rec(yN, y, deb, mil-1, opt)
    else
        return ajouter_rec(yN, y, mil+1, fin, opt)
    end
end

function ajouter!(yN::Vector{Solution}, 
                  y::Solution, 
                  opt::Optimisation=MAX)
    #! Affichage
    #println("Ajout de ", y)

    # Recherche de la position de y
    if length(yN) == 0
        insert!(yN, 1, y)
        #afficher(yN)
    else
        ind = ajouter_rec(yN, y, 1, length(yN), opt)
        if ind > 0 # Insertion
            insert!(yN, ind, y)
            #afficher(yN)

            # Elimination des points dominés
            verifier(yN, y, ind, opt)
            #println("Vérification : ")
            #afficher(yN)
        end
    end
end

## AFFICHER
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
