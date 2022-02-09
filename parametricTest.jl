################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Calcul paramétrique de la relaxation continue - Tests                #
################################################################################

include("dataStructures.jl")
include("functions.jl")
include("parametricMethod.jl")

using Test 

# ---------------------------------------------------------------------------- #
function testRatios()
	didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
	r1, r2 = ratios(didactic) 
	
	@testset "Ratio Tests" begin 
		@test r1 == [11//4, 1//2, 4//3, 5//2, 3//1, 1//2]
		@test r2 == [1//2, 7//4, 4//3, 1//1, 1//3, 3//2]
	end; 
end 

# ---------------------------------------------------------------------------- #
function testCriticalWeights() 
	didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
	r1, r2 = ratios(didactic) 
	weights, pairs = criticalWeights(didactic, r1, r2) 
	
	@testset "Critical Weights | Didactic" begin 
		@test pairs == [(1,4), (4,5), (1,5), (3,5), (1,3), (2,5), (1,2), (2,3), 
						(5,6), (1,6), (2,4), (3,4), (4,6), (3,6)] 
	end
end

# ---------------------------------------------------------------------------- #
function testTranspositions() 
	weights = [9//10, 85//100, 8//10, 8//10, 66//100, 65//100, 4//10, 4//10, 4//10]
	pairs = [(1,5), (2,4), (2,5), (3,6), (4,5), (1,2), (1,4), (5,6), (3,5)]
	@assert length(weights) == length(pairs) "There must be as many item pairs as there are critical weights"
	transpositions = transpositionPreprocessing(weights, pairs)
	
	@testset "Transposition Tests" begin 
		@test [t.λ for t in transpositions] == [9//10, 85//100, 8//10, 66//100, 65//100, 4//10]
		@test transpositions[1].pairs == [(1,5)] 
		@test transpositions[3].pairs == [(2,5), (3,6)] 
		@test transpositions[6].pairs == [(1,4), (5,6), (3,5)] 
	end;
end 

# ---------------------------------------------------------------------------- #
function testSwaps() 
	sequence = [3,1,2] 
	t = Transposition(1//2, [(1,2), (1,3), (2,3)]) 
	
	seq, pos = swaps3(sequence, t) 
	
	@testset "Swap Tests" begin
		@test seq == [2,1,3] 
		@test pos == sortperm(seq) 
	end;
end
