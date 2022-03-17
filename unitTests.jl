################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Unit tests                                                           #
################################################################################

include("lpRelaxation.jl")
include("martelloAndToth.jl")
include("branch-and-bound.jl")
using Test 

#prob = _MOMKP([11 2 2 8 10 9 1 9 16 ; 2 7 7 8 4 1 3 1 4], [4 4 4 6 4 3 2 3 6], [13])

# ----- INITIALISATION AND VARIABLE SETTING ---------------------------------- #
function testTranspositions() 
	weights = [9//10, 85//100, 8//10, 8//10, 66//100, 65//100, 4//10, 4//10, 4//10]
	pairs   = [(1,5), (2,4), (2,5), (3,6), (4,5), (1,2), (1,4), (5,6), (3,5)]

	@assert length(weights) == length(pairs) "There must be as many item pairs as there are critical weights"
	transpositions = transpositionPreprocessing(weights, pairs)
	
	@testset "Transposition Tests" begin 
		@test [t.λ for t in transpositions] == [9//10, 85//100, 8//10, 66//100, 65//100, 4//10]
		@test transpositions[1].pairs == [(1,5)] 
		@test transpositions[3].pairs == [(2,5), (3,6)] 
		@test transpositions[6].pairs == [(1,4), (5,6), (3,5)] 
	end;
end 

# Unit tests for setVariable
function testSetVariable()
	didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
	init     = initialisation(didactic, PARAMETRIC_LP)
	newInit  = setVariable(init, 4, PARAMETRIC_LP)

	@testset "SetVariable Tests" begin 
		@test newInit.seq == [5,1,3,2,6]
		@test newInit.pos == [2,4,3,3,1,5]
	end;
end

# ----- SOLUTIONS ------------------------------------------------------------ #
# Unit tests for reoptSolution
function testReoptSolution() 
	didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
	prev     = [5,1,4,3,2,6] 
	seq      = [5,3,4,1,2,6] 
	sol      = Solution{Float64}([1,0,0,1,1,0], [30,7], 0) 
	sol, s   = reoptSolution(didactic, seq, 2, 4, sol)
	
	@testset "Reoptimisation Tests" begin 
		@test sol.X  == [0,0,1,1//2,1,0]
		@test sol.z  == [22,11]
		@test sol.ω_ == 2
		@test s      == 3
	end;
	
	sol, s = reoptSolution(didactic, prev, 2, s, sol)
	
	@testset "Reoptimisation Tests 2" begin 
		@test sol.X  == [1,0,0,1,1,0]
		@test sol.z  == [30,7]
		@test sol.ω_ == 0
		@test s      == 4
	end;

	fname = "../instancesPG/set2/D4.DAT" 
	prob  = readInstanceMOMKPformatPG(false, fname)

	prev = [27, 15, 3, 6, 10, 28, 12, 13, 18, 23, 8, 2, 24, 16, 19, 26, 11, 9, 4, 
		20, 7, 17, 1, 14, 21, 5, 25, 22]
	seq  = [27, 15, 3, 6, 10, 28, 12, 13, 18, 23, 8, 2, 24, 16, 19, 26, 11, 20, 
		9, 4, 7, 17, 1, 14, 21, 5, 25, 22]
	sol  = Solution{Float64}(
		[0//1, 1//1, 1//1, 29//133, 0//1, 1//1, 0//1, 1//1, 1//1, 1//1, 1//1, 
		1//1, 1//1, 0//1, 1//1, 1//1, 0//1, 1//1, 1//1, 0//1, 0//1, 0//1, 1//1, 
		1//1, 0//1, 1//1, 1//1, 1//1], [135089//133, 147296//133], 29)
	
	sol, s = reoptSolution(prob, seq, 18, 19, sol)
	
	@testset "Reoptimisation Tests 3" begin 
		@test sol.z  == [135039//133, 155049//133]
		@test s      == 19 
		@test sol.ω_ == 73 
	end;	
end

# ----- MARTELLO AND TOTH ---------------------------------------------------- #
# Unit tests for chooseBound
function testChooseBound()
	# Didactic example 
	prob       = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
	r1, r2     = ratios(prob)
	weights, _ = criticalWeights(prob, r1, r2) 
	 
	upperBound = DualBoundSet{Float64}()
	
	@testset "chooseBound Tests" begin 
	
		# Case where λeq is smaller than the next critical weight   
		U0 = [30.0, 7.0]
		U1 = [23.0, 9.0]
		chooseBound!(upperBound, 1//1, weights[1], U0, U1)
		@test upperBound == [[30.0, 7.0]] 
		
		# Case where λeq is equal the next critical weight
		U0 = [30.0, 7.0]
		U1 = [20.0, 13.0]
		chooseBound!(upperBound, weights[3], weights[4], U0, U1)
		@test upperBound == [[30.0, 7.0]]
		
		# Case where λeq is in ]next, prev[
		# Cas λeq strictement compris entre les λ critiques précédent et suivant
		U0 = [30.0, 7.0]
		U1 = [83/4, 25/2]
		chooseBound!(upperBound, weights[4], weights[5], U0, U1)
		@test upperBound == [[83/4, 25/2], [30.0, 7.0]]
		
		
		U0 = [37/2, 55/4] 
		U1 = [25.0, 10.0] 
		chooseBound!(upperBound, weights[6], weights[7], U0, U1)
		@test upperBound == [[37/2, 55/4], [83/4, 25/2], [30.0, 7.0]]
	end;
end 

# ----- BRANCH-AND-BOUND ----------------------------------------------------- # 
# TODO : Unit tests for localNadirPoints and shiftedLocalNadirPoints

# Unit tests for isDominated 
function testIsDominated()

	UB = DualBoundSet{Float64}([[1.,8.], [4.,7.], [7.,4.], [9.,1.]],
		[Constraint(0//1, [1.,8.]), Constraint(7//15, [4.,7.]),
		 Constraint(3//5, [7.,4.]), Constraint(1//1, [9.,1.])])

	L1 = [Solution{Float64}(Rational{Int}[], [3.,9.], 0),
		  Solution{Float64}(Rational{Int}[], [6.,8.], 0), 
		  Solution{Float64}(Rational{Int}[], [9.,6.], 0),
		  Solution{Float64}(Rational{Int}[], [10.,4.], 0)]

	L2 = [Solution{Float64}(Rational{Int}[], [4.,9.], 0),
		  Solution{Float64}(Rational{Int}[], [8.,5.], 0),
		  Solution{Float64}(Rational{Int}[], [10.,3.], 0)]

	@testset "Dominance Tests" begin 
		@test isDominated(UB, L1)
		@test !isDominated(UB, L2)
	end;
end 

# ----- Main ----------------------------------------------------------------- #
# Main function that calls all the unit tests 
function mainTest()
	# Initialisation and variable setting 
	testTranspositions()
	testSetVariable()

	# Solutions 
	testReoptSolution() 

	# Martello and Toth 
	testChooseBound()

	# Branch-and-bound 
	testIsDominated()
end 

