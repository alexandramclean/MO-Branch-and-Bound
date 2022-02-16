################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Calcul paramétrique de la borne de Martello et Toth - Tests          #
################################################################################

include("dataStructures.jl")
include("functions.jl")
include("martelloAndToth.jl")

using Test 

# ---------------------------------------------------------------------------- #
function testChooseBound()
	# Exemple didactique
	prob = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
	
	r1, r2 = ratios(prob)
	weights, pairs = criticalWeights(prob, r1, r2) 
	 
	upperBound = Vector{Float64}[]
	
	@testset "ChooseBound Tests" begin 
		U0 = [30.0, 7.0]
		U1 = [23.0, 9.0]
		chooseBound!(upperBound, weights, U0, U1, 0)
		@test upperBound == [[30.0, 7.0]] 
		
		U0 = [30.0, 7.0]
		U1 = [21.5, 12.0]
		chooseBound!(upperBound, weights, U0, U1, 1)
		@test upperBound == [[30.0, 7.0]] 
		
		U0 = [30.0, 7.0]
		U1 = [83/4, 25/2]
		chooseBound!(upperBound, weights, U0, U1, 4)
		@test upperBound == [[83/4, 25/2], [30.0, 7.0]]
		
		
		U0 = [37/2, 55/4] 
		U1 = [25.0, 10.0] 
		chooseBound!(upperBound, weights, U0, U1, 6)
		@test upperBound == [[37/2, 55/4], [83/4, 25/2], [30.0, 7.0]]
	end;
	
end 
