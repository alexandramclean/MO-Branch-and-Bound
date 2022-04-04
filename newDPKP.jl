#using JuMP, GLPK

struct data2O1DKP
    p1::Vector{Int64} # profits objectif 1
    p2::Vector{Int64} # profits objectif 2
    w::Vector{Int64} # poids dimension 1
    Omega::Int64 # capacité dimension 1
    nbItems::Int64 # nombre d'objets
end

mutable struct ratios
    r1::Vector{Rational{Int64}} # ratio profit1/poids
    r2::Vector{Rational{Int64}} # ratio profit2/poids
end

mutable struct critique
    item1::Int64 # Objet 1 de la paire
    item2::Int64 # Objet 2 de la paire
    poids::Rational{Int64} # Poids critique de changement d'ordre des deux objets par rapport aux ratios pondéres
end

mutable struct objet
    indice::Int64 # Indice initial
    w::Int64 # poids
    p1::Int64 # profit 1
    p2::Int64 # profit 2
    r1::Rational{Int64} # ratio 1
    r2::Rational{Int64} # ratio 2
end

mutable struct extreme
    z1::Int64
    z2::Int64
    lambda::Rational{Int64}
end

struct solutionKP
    x::Vector{Int8}
    z1::Int64
    z2::Int64
end

#= Structure correspondant à un état basique (sans information de bornes...) pour une version basique de l'algorithme,
   nous n'avons ici qu'une solution partielle qu'on peut reconstruire par état =#
mutable struct etatBasic
	z1::Int64
    z2::Int64
	omega::Int64
	prec::Union{etatBasic,Nothing}
end


function parseKP(filename::String)
	f::IOStream = open(filename,"r")

    # Première ligne : taille du problème (nombre d'objets)
    s::String = readline(f)
    tab::Vector{Int64} = parse.(Int64,split(s," ",keepempty = false))
    @inbounds n::Int64 = tab[1]

	# Deuxième ligne : capacité du sac à dos
    s = readline(f)
	tab = parse.(Int64,split(s," ",keepempty = false))
	@inbounds Omega::Int64 = tab[1]

	# Troisième ligne : coefficients des coefficients de la première fonction objectif
    s = readline(f)
	p1::Vector{Int64} = parse.(Int64,split(s," ",keepempty = false))

    # Quatrième ligne : coefficients des coefficients de la seconde fonction objectif
    s = readline(f)
	p2::Vector{Int64} = parse.(Int64,split(s," ",keepempty = false))

	# Cinquième ligne : poids des objets
    s = readline(f)
	w::Vector{Int64} = parse.(Int64,split(s," ",keepempty = false))

    # End
    close(f)

    return data2O1DKP(p1,p2,w,Omega,n)
end

function parseKPbiVersMono(filename::String)
    f::IOStream = open(filename,"r")

    # Première ligne : taille du problème (nombre d'objets)
    s::String = readline(f)
    tab::Vector{Int64} = parse.(Int64,split(s," ",keepempty = false))
    @inbounds n::Int64 = tab[1]

    # Deuxième ligne : capacité du sac à dos
    s = readline(f)
	tab = parse.(Int64,split(s," ",keepempty = false))
	@inbounds Omega::Int64 = tab[1]

    # Troisième ligne ignorée
    s = readline(f)

    # Allocation mémoire pour les données
    p1::Vector{Int64} = Vector{Int64}(undef,n)
    p2::Vector{Int64} = Vector{Int64}(undef,n)
    w::Vector{Int64} = Vector{Int64}(undef,n)

    # Lignes suivantes
    for i in 1:n
        println(i)
        s = readline(f)
        tab = parse.(Int64,split(s,keepempty = false))
        p1[i] = tab[1]
        p2[i] = tab[2]
        w[i] = tab[3]
    end

    # End
    close(f)

    return data2O1DKP(p1,p2,w,Omega,n)
end

# Fonction traitant entièrement un nouvel état (teste sa dominance dans une colonne, et s'il n'est pas dominé filtre les états qu'il domine dans cette même colonne)
function traiteEtatBasic(z1part::Int64, z2part::Int64, wpart::Int64, etatPrec::etatBasic, colonneEtats::Vector{etatBasic}, taille::Int64,TableCap::Vector{Int64},tailleTableCap::Int64)
    if (taille == 0) # Il s'agit ici du premier état de la colonne #
        colonneEtats[taille + 1] = etatBasic(z1part,z2part,wpart,etatPrec)
        TableCap[1] = 1
        #println("insertion du premier etat")
        return 1,1
    else # Test de dominance du nouvel état nécessaire
        # Nous exploitons ici le fait que les états précédemment construits sont triés lexicographiquement pour chaque capacité
        domine::Bool = false
        i::Int64 = 1
        j::Int64 = 1
        while (i <= tailleTableCap) # Tant qu'on a des états dont les capacités n'interdisent pas qu'ils dominent le nouvel état
            j = TableCap[i]
            if (i == tailleTableCap)
                jmax = taille
            else
                jmax = TableCap[i+1] - 1
            end
            while (j <= jmax) && (colonneEtats[j].z1 >= z1part)
                if (colonneEtats[j].z2 >= z2part) # etat dominé
                    #println("($z1part,$z2part) -- $wpart est dominé par ($(colonneEtats[j].z1),$(colonneEtats[j].z2)) -- $(colonneEtats[j].omega)")
                    domine = true
                    j = jmax + 1 # Sortie de boucle
                else
                    #println("($z1part,$z2part) -- $wpart n'est pas dominé par ($(colonneEtats[j].z1),$(colonneEtats[j].z2)) -- $(colonneEtats[j].omega)")
                    j += 1
                end
            end
            if (domine)
                i = tailleTableCap + 1 # Sortie de boucle
            else
                i += 1
            end
        end
        if (!domine) # On crée l'état puisqu'il n'est pas dominé
            if (wpart != colonneEtats[taille].omega) # le remplissage augmente avec cet état => MAJ de TableCap
                TableCap[tailleTableCap + 1] = taille + 1
                tailleTableCap += 1
            end
            colonneEtats[taille + 1] = etatBasic(z1part,z2part,wpart,etatPrec)
            taille += 1
            #println("Ajout d'un état")
        end
        return taille,tailleTableCap
    end
end

function printColonne(C::Vector{etatBasic},taille::Int64)
    etat::etatBasic = C[1]
    for i in 1:taille
        etat = C[i]
        println("($(C[i].z1),$(C[i].z2)) -- $(C[i].omega)")
    end
    println("\n")
end

function DPBasic(d::data2O1DKP)
    # Par mesure de simplicité (pour éviter de taper d.X), quelques variables inutiles sont déclarées
    n::Int64 = d.nbItems
    Omega::Int64 = d.Omega
    p1::Vector{Int64} = d.p1
    p2::Vector{Int64} = d.p2
    w::Vector{Int64} = d.w

    # Simples variables de boucles
    i::Int64 = 0
    j::Int64 = 0

    #= Variable indiquant les poids des objets restants à considérer
       Utilisée uniquement dans le premier test de dominance (pas efficace mais pas cher...) =#
    omegarest::Int64 = 0
    for i in 1:n
        omegarest += w[i]
    end

    # Variables utilisées pour décrire le profit pour chaque fonction objectif et le poids d'une solution partielle dans l'algorithme de programmation dynamique
    z1part::Int64 = 0
    z2part::Int64 = 0
    wpart::Int64 = 0

    # L'ensemble des états est représenté par un tableau de tableau d'états
    # Un autre tableau est utilisé pour indiquer la taille utilisée (i.e. pas allouée) des tableaux d'états
    etats::Vector{Vector{etatBasic}} = Vector{Vector{etatBasic}}(undef,n + 1) # Tableau de tableaux d'états
    sizeEtats::Vector{Int64} = Vector{Int64}(undef,n + 1) # Taille des tableaux d'états

    # Quelques ajouts pratiques pour mieux parcourir les états et réduire un peu les coûts des tests de dominance
    TableCap::Vector{Int64} = Vector{Int64}(undef,Omega + 1) # Tableau qui sera utilisé pour indiquer les indices auxquels commence chaque nouvelle capacité résiduelle pour la colonne courante
    tailleTableCap::Int64 = 0 # taille du tableau ci-dessous
    omegaPrec::Int64 = 0 # remplissage pour les derniers états insérés (utile pour détecter les changements)

    # L'état initial est le seul de la première colonne
    etats[1] = Vector{etatBasic}(undef,1)
    etats[1][1] = etatBasic(0,0,0,nothing)
    sizeEtats[1] = 1

    sansobj::Int64 = 1
    avecobj::Int64 = 1 # Itérateurs pour quels états dans la colonne i (pendant la construction de la colonne i+1) on a construit l'état suivant sans objet ou avec objet (si possible)
    colonneprec::Vector{etatBasic} = etats[1]
    colonnesuiv::Vector{etatBasic} = etats[1] # Variable utilisées pour représenter la colonne précédante/suivante (utiles pour la lisibilité)
    #println("Affichage colonne 0")
    #printColonne(colonnesuiv,1)

    # Pour chaque colonne suivante (correspondant à l'inclusion de l'objet i-1), la dernière colonne est traitée séparément
    for i in 1:(n-1)
        #println("Remplissage colonne $i\n")
        # Initialisation de la colonne i + 1 des états
        etats[i + 1] = Vector{etatBasic}(undef,max(2*sizeEtats[i],Omega+1)) # Ici, il est choisi un dimensionnement dans le pire cas de la colonne courante pour éviter des réallocations
        sizeEtats[i + 1] = 0
        colonneprec = etats[i]
        colonnesuiv = etats[i+1]
        sansobj = 1
        avecobj = 1

        # On saute d'abord les états pour lesquels on met forcément une variable à 1 avec l'itérateur sansobj
        while (sansobj <= sizeEtats[i]) && (Omega - colonneprec[sansobj].omega >= omegarest)
            sansobj += 1
        end

        #println("sans obj init : $sansobj")
        while (sansobj <= sizeEtats[i])
            #println("etat potentiel sans obj : ($(colonneprec[sansobj].z1),$(colonneprec[sansobj].z2)) -- $(colonneprec[sansobj].omega)")
            #println("etat potentiel avec obj : ($(colonneprec[avecobj].z1 + p1[i]),$(colonneprec[avecobj].z2 + p2[i])) -- $(colonneprec[avecobj].omega + w[i])")
            # On cherche à chaque fois d'abord quel nouvel état donnerait la plus grande capacité (ou de manière équivalente le plus petit remplissage)
            if (colonneprec[sansobj].omega < colonneprec[avecobj].omega + w[i])
                #println("cas 1 : pas d'ajout car capacité inférieure")
                # On crée l'état pour lequel on n'ajoute pas d'objet s'il n'est pas dominé par un état existant
                sizeEtats[i + 1], tailleTableCap = traiteEtatBasic(colonneprec[sansobj].z1,colonneprec[sansobj].z2,colonneprec[sansobj].omega,colonneprec[sansobj],colonnesuiv,sizeEtats[i+1],TableCap,tailleTableCap)
                sansobj += 1
            elseif (colonneprec[sansobj].omega == colonneprec[avecobj].omega + w[i]) # Cas d'égalité
                if (colonneprec[sansobj].z1 > colonneprec[avecobj].z1 + p1[i])
                    #println("cas 2 : pas d'ajout car capacité égale mais z1 supérieur")
                    # On crée d'abord l'état pour lequel on n'ajoute pas d'objet s'il n'est pas dominé par un état existant
                    sizeEtats[i + 1], tailleTableCap = traiteEtatBasic(colonneprec[sansobj].z1,colonneprec[sansobj].z2,colonneprec[sansobj].omega,colonneprec[sansobj],colonnesuiv,sizeEtats[i+1],TableCap,tailleTableCap)
                    sansobj += 1
                elseif (colonneprec[sansobj].z1 < colonneprec[avecobj].z1 + p1[i])
                    #println("cas 3 : ajout car capacité égale mais z1 supérieur")
                    # On crée l'état pour lequel on ajoute un objet s'il n'est pas dominé par un état existant
                    sizeEtats[i + 1],tailleTableCap = traiteEtatBasic(colonneprec[avecobj].z1 + p1[i],colonneprec[avecobj].z2 + p2[i],colonneprec[avecobj].omega + w[i],colonneprec[avecobj],colonnesuiv,sizeEtats[i+1],TableCap,tailleTableCap)
                    avecobj += 1
                else
                    if (colonneprec[sansobj].z2 >= colonneprec[avecobj].z2 + p2[i])
                        #println("cas 4 : pas d'ajout car capacité égale, z1 égal mais z2 supérieur (double incrément)")
                        # On crée uniquement l'état pour lequel on n'ajoute pas d'objet s'il n'est pas dominé par un état existant
                        sizeEtats[i + 1], tailleTableCap = traiteEtatBasic(colonneprec[sansobj].z1,colonneprec[sansobj].z2,colonneprec[sansobj].omega,colonneprec[sansobj],colonnesuiv,sizeEtats[i+1],TableCap,tailleTableCap)
                        sansobj += 1
                        avecobj += 1
                    else
                        #println("cas 5 : ajout car capacité égale, z1 égal mais z2 supérieur (double incrément)")
                        # On crée uniquement l'état pour lequel on ajoute un objet s'il n'est pas dominé par un état existant
                        sizeEtats[i + 1],tailleTableCap = traiteEtatBasic(colonneprec[avecobj].z1 + p1[i],colonneprec[avecobj].z2 + p2[i],colonneprec[avecobj].omega + w[i],colonneprec[avecobj],colonnesuiv,sizeEtats[i+1],TableCap,tailleTableCap)
                        avecobj += 1
                        sansobj += 1
                    end
                end
            else # Ici, on a nécessairement que (colonneprec[sansobj].omega > colonneprec[avecobj].omega + w[i])
                #println("cas 6 : ajout car capacité inférieure")
                # On crée l'état pour lequel on ajoute un objet s'il n'est pas dominé par un état existant
                sizeEtats[i + 1],tailleTableCap = traiteEtatBasic(colonneprec[avecobj].z1 + p1[i],colonneprec[avecobj].z2 + p2[i],colonneprec[avecobj].omega + w[i],colonneprec[avecobj],colonnesuiv,sizeEtats[i+1],TableCap,tailleTableCap)
                avecobj += 1
            end
            #println("sansobj = $sansobj, avecobj = $avecobj")
            #readline(stdin)
        end

        # Il reste ensuite à avancer jusqu'à sizeEtats[i] pour l'itérateur avec objet
        # Il sera ici nécessaire de vérifier que cela mène bien à des solution admissibles
        while (avecobj <= sizeEtats[i]) && (colonneprec[avecobj].omega + w[i] <= Omega)
            #println("cas 7 : il n'y a plus que des cas d'ajouts à traiter")
            # On crée l'état pour lequel on ajoute un objet s'il n'est pas dominé par un état existant
            sizeEtats[i + 1],tailleTableCap = traiteEtatBasic(colonneprec[avecobj].z1 + p1[i],colonneprec[avecobj].z2 + p2[i],colonneprec[avecobj].omega + w[i],colonneprec[avecobj],colonnesuiv,sizeEtats[i+1],TableCap,tailleTableCap)
            avecobj += 1
            #println("avecobj = $avecobj")
        end
        omegarest -= w[i] # Mise à jour du poids des objets restants
        #println("Affichage colonne $i")
        #printColonne(colonnesuiv,sizeEtats[i+1])
        #if (i == 5)
        #    for i in 1:tailleTableCap
        #        println("$(TableCap[i])")
        #    end
        #end
        resize!(etats[i+1],sizeEtats[i+1]) # Ajustement de la taille de la colonne
    end

    # Traitement de la colonne finale
    etats[n + 1] = Vector{etatBasic}(undef,sizeEtats[n]) # Ici, il est choisi un dimensionnement dans le pire cas de la colonne courante pour éviter des réallocations
    tabTemp = Vector{etatBasic}(undef,sizeEtats[n])
    sizeEtats[n + 1] = 0
    colonneprec = etats[n]
    colonnesuiv = etats[n+1]

    # On fait d'abord une distinction entre les états pour lesquels il est possible d'ajouter un objet (et donc obligatoire)
    # et ceux pour lesquels un ajout est impossible
    i = 1
    while (i <= tailleTableCap) && (colonneprec[TableCap[i]].omega + w[n] <= Omega) i += 1 end
    sansobj = TableCap[i] # Indice du premier objet pour lequel on ne peut pas ajouter d'objet
    avecobj = sansobj - 1 # Indice du dernier objet pour lequel on peut ajouter d'objet

    # Un filtrage par dominance restera nécessaire pour les solutions obtenues
    # Pour cette dernière colonne, j'ai finalement choisi une solution simple (mais temporaire)
    # On instancie les états et on les trie lexicographiquement par rapport à (z1,z2) pour pouvoir gérer efficacement le test de dominance (malgré une recopie)

    j = 1
    for i in sansobj:sizeEtats[n]
        tabTemp[j] = etatBasic(colonneprec[i].z1,colonneprec[i].z2,colonneprec[i].omega,colonneprec[i])
        j += 1
    end
    for i in 1:avecobj
        tabTemp[j] = etatBasic(colonneprec[i].z1 + p1[n],colonneprec[i].z2 + p2[n],colonneprec[i].omega + w[n],colonneprec[i])
        j += 1
    end

    #println("Affichage colonne finale avant tri")
    #printColonne(tabTemp,sizeEtats[n])

    sort!(tabTemp,by = x -> 10000 * x.z1 + x.z2,rev = true)

    #println("Affichage colonne finale après tri")
    #printColonne(tabTemp,sizeEtats[n])

    colonnesuiv[1] = tabTemp[1]
    j = 1
    for i in 2:length(tabTemp)
        if (tabTemp[i].z2 > colonnesuiv[j].z2)
            colonnesuiv[j + 1] = tabTemp[i]
            j += 1
        end
    end
    sizeEtats[n+1] = j

    # Reconstruction des solutions en remontant le graphe jusqu'à l'état original
    X::Vector{Int8} = Vector{Int8}(undef,0)
    z1::Int64 = 0
    z2::Int64 = 0
    L::Vector{solutionKP} = Vector{solutionKP}(undef,sizeEtats[n+1])
    for i in 1:sizeEtats[n+1]
        courant::etatBasic = colonnesuiv[i]
        X = Vector{Int8}(undef,n)
        for i in (n + 1):-1:2
            if (courant.omega == (courant.prec).omega)
                X[i - 1] = Int8(0)
            else
                X[i - 1] = Int8(1)
            end
            courant = courant.prec
        end
        z1 = 0
        z2 = 0
        for i in 1:n
            z1 += p1[i]*X[i]
            z2 += p2[i]*X[i]
        end
        L[i] = solutionKP(X,z1,z2)
    end

    # Retour du tableau de solutions obtenues
    return L
end

function epsilonC(d::data2O1DKP)
    m = Model(GLPK.Optimizer)
    n::Int64 = d.nbItems
    @variable(m,x[1:n],binary = true)
    @constraint(m,dim1,sum(d.w[i]*x[i] for i in 1:n) <= d.Omega)
    @objective(m,Max,sum(d.p1[i]*x[i] for i in 1:n))

    # Première résolution
    optimize!(m)
    status = termination_status(m)
    X::Vector{Int64} = round.(Int64,value.(x))
    z1::Int64 = 0
    z2::Int64 = 0
    for i in 1:n
        z1 += d.p1[i]*X[i]
        z2 += d.p2[i]*X[i]
    end
    Ltemp::Vector{solutionKP} = Vector{solutionKP}(undef,0)
    push!(Ltemp,solutionKP(X,z1,z2))

    # Résolutions suivantes

    @variable(m, const_term)
    @constraint(m,obj2,sum(d.p2[i]*x[i] for i in 1:n) >= const_term + 1)
    fix(const_term,z2)
    optimize!(m)
    status = termination_status(m)

    while (status != MOI.INFEASIBLE)
        X = round.(Int64,value.(x))
        z1 = 0
        z2 = 0
        for i in 1:n
            z1 += d.p1[i]*X[i]
            z2 += d.p2[i]*X[i]
        end
        push!(Ltemp,solutionKP(X,z1,z2))
        fix(const_term,z2)
        optimize!(m)
        status = termination_status(m)
    end

    L::Vector{solutionKP} = Vector{solutionKP}(undef,0)
    for i in 1:(length(Ltemp) - 1)
        if (Ltemp[i].z1 != Ltemp[i+1].z1)
            push!(L,Ltemp[i])
        end
    end
    push!(L,Ltemp[end])

    return L
end

#=
# Fonction de résolution récursive pour la résolution dichotomique tenant compte du risque d'avoir des solutions initiales dominées
function dichoRec(d::data2O1DKP,tabSol::Vector{solutionKP},yr::Tuple{Int64,Int64},yl::Tuple{Int64,Int64},y1::Tuple{Int64,Int64},y2::Tuple{Int64,Int64},dom::Vector{Bool})
	# yl solution à "gauche" et yr solution à "droite"
	lambda::Tuple{Int64,Int64} = (yl[2]-yr[2],yr[1]-yl[1]) # coefficients de la normale (positive) à la pente
	xetoile::Vector{Int8},zetoile::Tuple{Int64,Int64} = getOptP(m,d,lambda[1],lambda[2]) # nouvelle solution
	comparaison::Int64 = LinearAlgebra.dot(lambda,zetoile) - LinearAlgebra.dot(lambda,yr)

	# Si la nouvelle solution est au-dessus de la droite, on la garde
	if comparaison > 0
        push!(tabSol,solutionInit(xetoile,zetoile))
        if (zetoile[1] == y1[1]) # Si on obtient une solution optimale pour z1, elle domine nécessairement la solution initiale car comparaison > 0
            dom[1] = true # la solution initiale est indiquée comme étant dominée pour une pas être ajoutée ultérieurement
		    tabSol = dichoRec2(m,d,tabSol,zetoile,yl,y1,y2,dom) # et on n'a qu'un seul appel récursif
        elseif (zetoile[2] == y2[2]) # Si on obtient une solution optimale pour z2, elle domine nécessairement la solution initiale car comparaison > 0
            dom[2] = true # la solution initiale est indiquée comme étant dominée pour une pas être ajoutée ultérieurement
            tabSol = dichoRec2(m,d,tabSol,yr,zetoile,y1,y2,dom) # et on n'a qu'un seul appel récursif
        else # Cas classique
            tabSol = dichoRec2(m,d,tabSol,zetoile,yl,y1,y2,dom)
            tabSol = dichoRec2(m,d,tabSol,yr,zetoile,y1,y2,dom)
        end
	# Si l'on a les mêmes valeurs mais des solutions différentes
	elseif comparaison == 0 && zetoile != yr && zetoile != yl
        push!(tabSol,solutionInit(xetoile,zetoile))
	# Si l'on a fini -> path relinking ?
	#else
		#pathrelinking(yl,yr) ?
	end
	return tabSol
end


# Méthode dichotomique ne démarrant pas nécessairement par une résolution lexicographique
# Avantage : si on a par chance une solution lexicographique, on fait une résolution de moins
# Inconvénient : C'est un peu plus compliqué!
function dichoStart(d::donnees2O2DKP)
    # Déclaration d'un tableau de solutions
    tabSol::Vector{solutionInit} = Vector{solutionKP}(undef,0)

    # Résolution sur la première fonction objectif
    prob::Vector{item} = Vector{item}(undef,5)
    maxZ::Int64 = 0
    for i in 1:d.nbItems
        prob[i] = item(d.p1[i],d.omega[i],0,i-1)
        maxZ += d.p1[i]
    end
    z1::Int64 = Int64(comboJulia(prob,d.nbItems,d.Omega,maxZ))
    X::Vector{Int8} = Vector{Int8}(undef,d.nbItems)
    for i in 1:d.nbItems
        x[prob[i].i + 1] = prob[i].x
    end
    z2::Int64 =
    println("Valeur optimale = ",val)
    println("x = ",x)
    (x1::Vector{Int8},y1::Tuple{Int64,Int64}) = getOptP(m,d,1,0)

    # Résolution sur la seconde fonction objectif
    (x2::Vector{Int8},y2::Tuple{Int64,Int64}) = getOptP(m,d,0,1)

    # Lancement de la méthode dichotomique
    dom::Vector{Bool} = [false,false]
    dichoRec2(m,d,tabSol,y1,y2,y1,y2,dom)

    # Insertion des solutions initiales si elles ne sont pas dominées
    if (dom[1] == false)
        push!(tabSol,solutionInit(x1,y1))
    end
    if (dom[2] == false)
        push!(tabSol,solutionInit(x2,y2))
    end

    # Tri final (utile avant d'appliquer le path-relinking, ou une méthode en deux phases)
    sort!(tabSol, by = sol -> sol.z[1])

    # Retour final
    return tabSol
end
=#

#=
# Fonction de calculs des ratios (profit1/poids, profit2/poids) pour chaque objet
function calculRatios(d::donnees2O1DKP)
    r1::Vector{Rational{Int64}} = Vector{Rational{Int64}}(undef,d.nbItems)
    r2::Vector{Rational{Int64}} = Vector{Rational{Int64}}(undef,d.nbItems)
    for i in 1:d.nbItems
        r1[i] = d.p1[i]//d.w[i]
        r2[i] = d.p2[i]//d.w[i]
    end
    return ratios(r1,r2)
end

# Fonction calculant les poids critiques associés à chaque paire d'objets
# En ne conservant que ceux compris dans l'intervalle ]0,1[
# Puis en les triant par ordre décroissant
function calculPoidsCritiques(r::ratios)
    longueur::Int64 = length(r.r1)
    poids::Rational{Int64} = 0
    den::Rational{Int64} = 0
    crit::Vector{critique} = Vector{critique}(undef,0)
    for i in 1:(longueur - 1)
        for j in i+1:longueur
            den = r.r1[i] - r.r2[i] - r.r1[j] + r.r2[j]
            if (den != 0) # dénominateur nul => pas de poids critique
                poids = (r.r2[j] - r.r2[i]) // den
                if ((poids > 0) && (poids < 1)) # Seuls les changements dans l'intervalle ]0,1[ nous intéressent
                    push!(crit,critique(i,j,poids))
                end
            end
        end
    end
    sort!(crit, by = x->x.poids, rev = true)
    return crit
end

function RClambda1(d::donnees2O1DKP,r::ratios)
    nbItems::Int64 = d.nbItems
    Omega::Int64 = d.Omega
    omega::Int64 = 0
    objets::Vector{objet} = Vector{objet}(undef,d.nbItems)
    sequence::Vector{Int64} = Vector{Int64}(undef,d.nbItems)
    i::Int64 = 0
    s::Int64 = 0
    dantzig1::Int64 = 0
    dantzig2::Int64 = 0

    # Tri lexicographique des objets par ratios décroissants
    for i in 1:d.nbItems
        objets[i] = objet(i,d.w[i],d.p1[i],d.p2[i],r.r1[i],r.r2[i])
    end
    sort!(objets,by = x -> 1000*x.r1 + x.r2,rev = true)

    for i in 1:d.nbItems
        sequence[i] = objets[i].indice
    end

    # Construction de la solution dantzig
    i = 1
    omega += objets[i].w
    dantzig1 += objets[i].p1
    dantzig2 += objets[i].p2
    while (i <= d.nbItems) && (omega <= Omega) # Prise en compte des cas triviaux
        i += 1
        omega += objets[i].w
        dantzig1 += objets[i].p1
        dantzig2 += objets[i].p2
    end # Remplissage jusqu'à débordement du sac
    if (omega <= Omega)
        return objets,dantzig1,dantzig2,omega,s,true,true # Cas trivial d'une solution contenant tous les objets => Solution idéale du problème bi-objectif
    end
    s = i # s indique l'indice du break item
    dantzig1 = dantzig1 - objets[s].p1
    dantzig2 = dantzig2 - objets[s].p2
    omega = omega - objets[s].w
    if (omega == Omega)
        return objets,dantzig1,dantzig2,omega,s,true,false # Cas trivial où on obtient une solution admissible optimale => Nécessité de continuer la résolution
    else
        return objets,dantzig1, dantzig2,omega,s,false,false # Solution initiale contenant un fragment du break item
    end
end

function RCbiKP(d::donnees2O1DKP)
    # Calculs des ratios
    r::ratios = calculRatios(d)

    # Calculs des poids critiques
    crit::Vector{critique} = calculPoidsCritiques(r)

    # Calcul de la solution initiale
    objets::Vector{objet},dantzig1::Int64,dantzig2::Int64,omega::Int64,s::Int64,feasible::Bool,ideal::Bool = RClambda1(d,r)

    # Initialisation du tableau permettant de retrouver les objets
    indSequence::Vector{Int64} = Vector{Int64}(undef,d.nbItems)
    j::Int64 = 1
    for o in objets
        indSequence[o.indice] = j
        j += 1
    end

    # Initialisation du tableau de solutions admissible obtenues
    tabSol::Vector{solutionKP} = Vector{solutionKP}(undef,0)
    x::Vector{Int8} = Vector{Int8}(undef,d.nbItems)
    for i in 1:d.nbItems
        if (indSequence[i] < s)
            x[i] = 1
        else
            x[i] = 0
        end
    end

    if feasible == true
        push!(tabSol,solutionKP(x,dantzig1,dantzig2))


    # Lancement de la boucle principale
    if (ideal == true)
        # Créer l'unique solution obtenu en spécifiant qu'elle minimise les deux objectifs simultanément
    else
        for c in crit
            # Détermination d'un score pour savoir dans quel cas nous sommes
            score = 0
            ind1 = indSequence[c.item1]
            ind2 = indSequence[c.item2]
            if (ind1 < s) # item1 est dans le sac
                score += 2
            elseif (ind1 == s) && (!feasible) # item1 est le break item
                score += 1
            end
            if (ind2 < s) # item1 est dans le sac
                score += 2
            elseif (ind2 == s) && (!feasible) # item1 est le break item
                score += 1
            end

            # Swap à effectuer dans tous les cas
            tempo = objets[ind1]
            objets[ind1] = objets[ind2]
            objets[ind2] = tempo
            temp = ind1
            indSequence[c.item1] = ind2
            indSequence[c.item2] = temp

            if (score == 1) # Échange d'un objet pas dans le sac et du break item
                # On repart de la solution dantzig précédente
                i = s
                while (omega <= Omega)
                    i += 1
                    omega += objets[i].w
                    dantzig1 += objets[i].p1
                    dantzig2 += objets[i].p2
                end # Remplissage jusqu'à débordement du sac
                s = i # s indique l'indice du break item si la solution n'est pas admissible
                dantzig1 = dantzig1 - objets[s].p1
                dantzig2 = dantzig2 - objets[s].p2
                omega = omega - objets[s].w
                if (omega == Omega) # Test d'admissibilité
                    feasible = true
                else
                    feasible = false
                end
                # Solution à compléter ensuite
            elseif (score == 2) || (score = 3)
                # On repart de la solution dantzig à laquelle on enlève le dernier objet
                dantzig1 = dantzig1 - objets[s].p1
                dantzig2 = dantzig2 - objets[s].p2
                omega = omega - objets[s].w
                i = s-1
                while (omega <= Omega)
                    i += 1
                    omega += objets[i].w
                    dantzig1 += objets[i].p1
                    dantzig2 += objets[i].p2
                end # Remplissage jusqu'à débordement du sac
                s = i # s indique l'indice du break item si la solution n'est pas admissible
                dantzig1 = dantzig1 - objets[s].p1
                dantzig2 = dantzig2 - objets[s].p2
                omega = omega - objets[s].w
                if (omega == Omega) # Test d'admissibilité
                    feasible = true
                else
                    feasible = false
                end
            end
        end
    end
end

function essai()
    d::donnees2O1DKP = donnees2O1DKP([11,2,8,10,9,1],[2,7,8,4,1,3],[4,4,6,4,3,2],11,6)
    RCbiKP(d)
end
=#
