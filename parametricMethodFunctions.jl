################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Fonctions auxiliaires                                                #
################################################################################

include("functions.jl")
using TimerOutputs

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

