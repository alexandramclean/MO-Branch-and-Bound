################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Calcul paramétrique - Tests                                          #
################################################################################

include("lpRelaxation.jl")
include("martelloAndToth.jl")
using Test 

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
	sol, s, residualCapacity = reoptSolution(didactic, seq, 2, 4, sol, 0)
	
	@testset "Reoptimisation Tests" begin 
		@test sol.X == [0,0,1,1//2,1,0]
		@test sol.z == [22,11]
		@test residualCapacity == 2
		@test s == 3
	end;
	
	sol, s, residualCapacity = reoptSolution(didactic, prev, 2, s, sol, residualCapacity)
	
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
	sol, s, residualCapacity = reoptSolution(prob, seq, 18, 19, sol, 29)
	
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
	
# ---------------------------------------------------------------------------- #
function testChooseBound()
	# Exemple didactique
	prob = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
	
	r1, r2 = ratios(prob)
	weights, pairs = criticalWeights(prob, r1, r2) 
	 
	upperBound = Vector{Float64}[]
	constraints = Constraint[]
	
	@testset "chooseBound Tests" begin 
	
		# Cas λeq plus petit que le λ critique suivant 
		U0 = [30.0, 7.0]
		U1 = [23.0, 9.0]
		chooseBound!(upperBound, constraints,1//1, weights[1], U0, U1)
		@test upperBound == [[30.0, 7.0]] 
		
		# Cas λeq égal au λ critique suivant
		U0 = [30.0, 7.0]
		U1 = [20.0, 13.0]
		chooseBound!(upperBound, constraints, weights[3], weights[4], U0, U1)
		@test upperBound == [[30.0, 7.0]]
		
		# Cas λeq strictement compris entre les λ critiques précédent et suivant
		U0 = [30.0, 7.0]
		U1 = [83/4, 25/2]
		chooseBound!(upperBound, constraints, weights[4], weights[5], U0, U1)
		@test upperBound == [[83/4, 25/2], [30.0, 7.0]]
		
		
		U0 = [37/2, 55/4] 
		U1 = [25.0, 10.0] 
		chooseBound!(upperBound, constraints, weights[6], weights[7], U0, U1)
		@test upperBound == [[37/2, 55/4], [83/4, 25/2], [30.0, 7.0]]
	end;
end 



