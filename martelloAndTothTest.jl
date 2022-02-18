################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Calcul paramétrique de la borne de Martello et Toth - Tests          #
################################################################################

include("dataStructures.jl")
include("functions.jl")
include("martelloAndToth.jl")
include("parserMomkpPG.jl")
include("parserMomkpZL.jl")

using Test 

# ---------------------------------------------------------------------------- #
function testChooseBound()
	# Exemple didactique
	prob = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
	
	r1, r2 = ratios(prob)
	weights, pairs = criticalWeights(prob, r1, r2) 
	 
	upperBound = Vector{Float64}[]
	constraints = Constraint[]
	
	@testset "ChooseBound Tests" begin 
	
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

# ---------------------------------------------------------------------------- #
function testCasEgalite()
	# Exemple pathologique
	prob = _MOMKP([11 2 2 8 10 9 1 9 16 ; 2 7 7 8 4 1 3 1 4], [4 4 4 6 4 3 2 3 6], [13])
	
	#=r1, r2 = ratios(prob)
	weights, pairs = criticalWeights(prob, r1, r2)
	transpositions = transpositionPreprocessing(weights, pairs)
	
	seq = sortperm(1000000*r1 + r2, rev=true) # Item sequence 
	pos = sortperm(seq)          			  # Item positions
	println("seq : ", seq)
	
	sol, s, ω_ = dantzigSolution(prob, seq)
	println("X = ", sol.X, " z = ", sol.z)
	println("s = ", s, " ω_ = ", ω_)
	
	for t in transpositions
		println("λ = ", t.λ, " pairs = ", t.pairs)
	end =#
	
	upperBound, constraints = martelloAndToth(prob)

end

# ---------------------------------------------------------------------------- #
function testMartelloAndToth(dir::String)

	# Exemple didactique
	println("Exemple didactique")
	prob = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
	@time upperBound, constraints = martelloAndToth(prob)
	
	# Test sur toutes les instances du répertoire dir
	files = readdir(dir)
	for f in files 
	
		println("\n", f)
		
		fname = dir*f
		if fname[length(fname)-3:length(fname)] == ".DAT"
			prob = readInstanceMOMKPformatPG(false, fname)
		else
			prob = readInstanceMOMKPformatZL(false, fname)
		end
		
		@time upperBound, constraints = martelloAndToth(prob)
	end 

end
