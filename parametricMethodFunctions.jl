################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Fonctions auxiliaires                                                #
################################################################################

include("functions.jl")

# ----- RATIOS AND CRITICAL WEIGHTS ------------------------------------------ #
# Computes the critical weights
function criticalWeights(prob::_MOMKP,
                         r1::Vector{Rational{Int}},
                         r2::Vector{Rational{Int}})

    n       = size(prob.P)[2]
    weights = Rational{Int}[]
    pairs   = Tuple{Int,Int}[]

    nbTransp = 0

    # Computes the critical weight for each pair of items (i,j)
    for i in 1:n
        for j in i+1:n

            if !(r1[i] == r1[j] || r2[i] == r2[j]) 

                λ = (r2[j] - r2[i])//(r1[i]-r2[i]-r1[j]+r2[j])

                if λ > 0 && λ < 1
                    nbTransp += 1 
                    push!(weights, λ)
                    push!(pairs, (i,j))
                end
            end
        end
    end

    println("\tNumber of transpositions : ", nbTransp, " / ", Int(n*(n-1)/2))

    # Sorts the critical weights and associated item pairs in decreasing order
    perm = sortperm(weights, rev=true)
    return weights[perm], pairs[perm]
end

# Returns a list containing all the distinct λ values in decreasing order and
# the associated item pair(s)
function transpositionPreprocessing(weights::Vector{Rational{Int}},
                                    pairs::Vector{Tuple{Int,Int}})

    transpositions = Transposition[]

    iter = 1
    while iter <= length(weights)

        # There are multiple identical critical weights λ
        if iter < length(weights) && weights[iter] == weights[iter+1]

            transp = Transposition(weights[iter])

            while iter < length(weights) && weights[iter] == weights[iter+1]
            push!(transp.pairs, pairs[iter])
            iter += 1
            end

            # The last occurence of λ
            push!(transp.pairs, pairs[iter])
            iter += 1

            push!(transpositions, transp)

        else
            push!(transpositions, Transposition(weights[iter], [pairs[iter]]))
            iter += 1
        end
    end

    return transpositions
end

# ----- SEQUENCE AND POSITIONS ----------------------------------------------- #
# Updates the positions after reversing the subsequence between start and finish
function updatePositions!(seq::Vector{Int}, 
                         pos::Vector{Int}, 
                         start::Int, 
                         finish::Int) 

    if finish - start >= 1
        # Swaps the positions of the items at the edge of the subsequence 
        tmp              = pos[seq[start]] 
        pos[seq[start]]  = pos[seq[finish]]
        pos[seq[finish]] = tmp 

        # Appel récursif sur la sous-séquence privée des extrémités
        updatePositions!(seq, pos, start+1, finish-1)
    end
end

# ----- INITIALISATION ------------------------------------------------------- #
function initialisation(prob::_MOMKP)

	# Ratios and critical weights
	r1, r2 = ratios(prob)
	weights, transpositions = criticalWeights(prob, r1, r2)
	
	# Identical critical weights are grouped together 
	transpositions = transpositionPreprocessing(weights, pairs)
	
	# Sorts the ratios in lexicographically decreasing order according to (r1, r2) 
	seq = sortperm(1000000*r1 + r2, rev=true) # Item sequence 
	pos = sortperm(seq)          			  # Item positions
	
	return transpositions, seq, pos
end