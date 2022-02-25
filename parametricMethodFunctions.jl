################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Fonctions auxiliaires                                                #
################################################################################

include("functions.jl")
using TimerOutputs

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
# Identification of the modified subsequences
function identifySubsequences(positions)
    subsequences = Tuple{Int,Int}[]
	start = positions[1][1] ; finish = positions[1][2]

	for p in positions[2:end]

		if p[1] > finish # Start of a new distinct subsequence
			push!(subsequences, (start, finish))
			start = p[1] ; finish = p[2]
		else
			finish = p[2]
		end
	end
	push!(subsequences, (start, finish))
    return subsequences
end

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

    @timeit to "Initialisation" begin 
        # Ratios and critical weights
        @timeit to "Ratios" r1, r2 = ratios(prob)
        @timeit to "Critical weights" weights, pairs = criticalWeights(prob, r1, r2)
        
        # Identical critical weights are grouped together 
        @timeit to "Transpositions" transpositions = transpositionPreprocessing(weights, pairs)
        
        # Sort the ratios in lexicographically decreasing order according to (r1,r2) 
        @timeit to "Sequence" seq = sortperm(1000000*r1 + r2, rev=true) # Item sequence 
        @timeit to "Positions" pos = sortperm(seq)          			# Item positions
	end
	return transpositions, seq, pos
end

# ----- SETTING VARIABLES ---------------------------------------------------- #
# Remove a variable from the sequence and transpositions
function setVariable(transpositions::Vector{Transposition},
                     seq::Vector{Int},
                     pos::Vector{Int},
                     var::Int)
       
    newTranspositions = Transposition[]
    newSeq            = Vector{Int}(undef,length(seq)-1)
    newPos            = copy(pos)

    # The variable is removed from the set of transpositions
    for t in transpositions

        if length(t.pairs) > 1	

            swaps = Tuple{Int,Int}[]
            for pair in t.pairs 
                if !(var in pair) 
                    push!(swaps, pair)
                end
            end

            if length(swaps) > 0 
                push!(newTranspositions, Transposition(t.λ, swaps))
            end
        else

            if !(var in t.pairs[1])
                push!(newTranspositions, Transposition(t.λ, t.pairs))
            end
        end
    end

    # The variable is removed from the sequence
    inser = 1
    for i in 1:length(seq)
        if seq[i] != var
            newSeq[inser] = seq[i] 
            inser += 1 
        end
    end

    # The positions of items after var in the sequence are diminished by 1
    for p in pos[var]+1:length(pos)
        newPos[seq[p]] = pos[seq[p]] - 1
    end

    return newTranspositions, newSeq, newPos
end