################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Calcul paramétrique de la relaxation continue - Tests                #
################################################################################

include("lpRelaxation.jl")

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
function testReoptSolution() 
	didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
	prev = [5,1,4,3,2,6] 
	seq  = [5,3,4,1,2,6] 
	sol  = Solution([1,0,0,1,1,0], [30,7]) 
	sol, s, residualCapacity = reoptSolution(didactic, prev, seq, 2, 4, sol, 0)
	
	@testset "Reoptimisation Tests" begin 
		@test sol.X == [0,0,1,1//2,1,0]
		@test sol.z == [22,11]
		@test residualCapacity == 2
		@test s == 3
	end;
	
	sol, s, residualCapacity = reoptSolution(didactic, seq, prev, 2, s, sol, residualCapacity)
	
	@testset "Reoptimisation Tests 2" begin 
		@test sol.X == [1,0,0,1,1,0]
		@test sol.z == [30,7]
		@test residualCapacity == 0
		@test s == 4
	end;
	
	prev = [27, 15, 3, 6, 10, 28, 12, 13, 18, 23, 8, 2, 24, 16, 19, 26, 11, 9, 4, 20, 7, 17, 1, 14, 21, 5, 25, 22]
	seq = [27, 15, 3, 6, 10, 28, 12, 13, 18, 23, 8, 2, 24, 16, 19, 26, 11, 20, 9, 4, 7, 17, 1, 14, 21, 5, 25, 22]
	sol = Solution([0//1, 1//1, 1//1, 29//133, 0//1, 1//1, 0//1, 1//1, 1//1, 1//1, 1//1, 1//1, 1//1, 0//1, 1//1, 1//1, 0//1, 1//1, 1//1, 0//1, 0//1, 0//1, 1//1, 1//1, 0//1, 1//1, 1//1, 1//1], [135089//133, 147296//133])
	
	fname = "../instancesPG/set2/D4.DAT" 
	prob  = readInstanceMOMKPformatPG(false, fname)
	sol, s, residualCapacity = reoptSolution(prob, prev, seq, 18, 19, sol, 29)
	
	@testset "Reoptimisation Tests 3" begin 
		@test sol.z == [135039//133, 155049//133]
		@test s == 19 
		@test residualCapacity == 73 
	end;	
end

function testSetVariable()
	didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
	transpositions, seq, pos = initialisation(didactic)
	newTranspositions, newSeq, newPos = setVariable(transpositions, seq, pos, 4)
	println(newSeq)
	println(newPos)
end
	



