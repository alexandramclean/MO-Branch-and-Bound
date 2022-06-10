################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Tests et comparaison de différentes méthodes pour la relaxation      #
################################################################################

include("lpRelaxation.jl")
include("martelloAndToth.jl")
include("dichotomicMethod.jl")
include("simplexAlgorithm.jl")
include("branch-and-bound.jl")
include("efficientSupportedSolutions.jl")
include("vOptMomkp.jl")
include("parsers.jl")
include("displayGraphic.jl")

using TimerOutputs
const to = TimerOutput()

# ----- DICHOTOMIC METHOD ---------------------------------------------------- #
# Compares the parametric method and the dichotomic method for computing the 
# LP relaxation of instance prob
# If graphic=true, produces a figure showing both obtained upper bound sets 
function compareParametric_Dichotomic(name, prob, graphic=false)

	@timeit to "\nParametric v. Dichotomic" begin 
		println("Parametric method")
		@timeit to "Parametric method" begin 
			@timeit to "Initialisation" init = initialisation(prob, PARAMETRIC_LP)
			@timeit to "Setvar" setvar = initialSetvar(prob, init, PARAMETRIC_LP)
			@timeit to "Relaxation" UBparam, _ = 
				parametricMethod(prob, init, setvar)
		end 

		println("Dichotomic method")
		@timeit to "Dichotomic method" begin 
			@timeit to "Initialisation" initDicho = initialisation(prob, DICHOTOMIC)
			@timeit to "Relaxation" UBdicho = 
				dichotomicMethod(prob, initDicho)
		end 
	end 

	if graphic 
   		# Setup
    	figure("Parametric method and dichotomic method | "*name,figsize=(6.5,5))
    	xlabel(L"z^1(x)")
    	ylabel(L"z^2(x)")
    	PyPlot.title("LP Relaxation | "*name)

		# Show the points computed by the parametric method 
		y_PN11 = [y[1] for y in UBparam.points] 
		y_PN12 = [y[2] for y in UBparam.points]
		scatter(y_PN11, y_PN12, color="green", marker="+", label = "parametric")
		plot(y_PN11, y_PN12, color="green", linewidth=0.75, marker="+",
			markersize=1.0, linestyle=":")

		# Show the points computed by the dichotomic method 
		y_PN21 = sort([y.point[1] for y in UBdicho]) 
		y_PN22 = sort([y.point[2] for y in UBdicho],rev=true)
		scatter(y_PN21, y_PN22, color="red", marker="+", label = "dichotomic")
		plot(y_PN21, y_PN22, color="red", linewidth=0.75, marker="+",
			markersize=1.0, linestyle=":")

    	legend(bbox_to_anchor=[1,1], loc=0, borderaxespad=0, fontsize = "x-small")
    end 
end

# Calls compareParametric_Dichotomic on all instances in directory dir 
function testInstances(dir::String, graphic=false) 

	println("Exemple didactique")
	didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])

	init   = initialisation(didactic, PARAMETRIC_LP)
	setvar = initialSetvar(didactic, init, PARAMETRIC_LP)
	_ = parametricMethod(didactic, init, setvar)

	initDicho = initialisation(didactic, DICHOTOMIC)
	_ = dichotomicMethod(didactic, initDicho)
	
	files = readdir(dir*"dat/") 
	for fname in files 
		println("\n", fname) 
		if fname[length(fname)-3:length(fname)] == ".DAT"
			prob = readInstanceMOMKPformatPG(false, dir*"dat/"*fname)
		elseif fname[length(fname)-3:length(fname)] == ".dat"
			prob = readInstanceKP(dir*"dat/"*fname)
		else
			prob = readInstanceMOMKPformatZL(false, dir*"dat/"*fname)
		end
	
		compareParametric_Dichotomic(fname, prob, graphic)
	end	
end 


# ----- MARTELLO AND TOTH ---------------------------------------------------- #
# Compare CPU times for the LP relaxation and the Martello and Toth upper bound
function compareLP_MT(prob::_MOMKP)

	@timeit to "\nLP v. MT" begin 
		# Initialisation
		init = initialisation(prob, PARAMETRIC_LP)
		setvar = initialSetvar(prob, init, PARAMETRIC_LP)

		println("LP Relaxation")
		@timeit to "LP Relaxation" UBparam = 
			parametricMethod(prob, init, setvar)
		
		println("Martello and Toth")
		@timeit to "Martello and Toth" UB = 
			martelloAndToth(prob, init, Solution{Float64}(prob))
	end 
end 

# Calls compareLP_MT on all instances in dir
function testInstancesMT(dir::String)

	println("Exemple didactique")
	didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])

	init   = initialisation(didactic, PARAMETRIC_LP)
	setvar = initialSetvar(didactic, init, PARAMETRIC_LP)

	_ = parametricMethod(didactic, init, setvar) 
	_ = martelloAndToth(didactic, init, Solution{Float64}(didactic))
	
	files = readdir(dir)
	for fname in files
		println("\n", fname)
		
		if fname[length(fname)-3:length(fname)] == ".DAT"
			prob = readInstanceMOMKPformatPG(false, dir*fname)
		else
			prob = readInstanceMOMKPformatZL(false, dir*fname)
		end
		println("n = ", size(prob.P)[2])
		compareLP_MT(prob)
	end
end

# ----- SETTING VARIABLES ---------------------------------------------------- #
# TODO : Compare CPU times for the initialisation and setVariable functions 


# ----- SIMPLEX ALGORITHM ---------------------------------------------------- #
# Compares the parametric method and the simplex algorithm for computing the 
# LP relaxation of instance prob
# If graphic=true, produces a figure showing both obtained upper bound sets 
function compareParametric_Simplex(prob::_MOMKP, fname::String, graphic=false)

	#newProb = groupEquivalentItems(prob)
	@timeit to "\nParametric v. Simplex" begin 
		println("Parametric method")
		@timeit to "Parametric method" begin 
			@timeit to "Initialisation" init   = initialisation(prob, PARAMETRIC_LP)
			@timeit to "Initial setvar" setvar = initialSetvar(prob, init, PARAMETRIC_LP)
			@timeit to "LP Relaxation" UBparam, _ = parametricMethod(prob, init, setvar)
		end 

		println("Simplex algorithm")
		@timeit to "Simplex algorithm" begin 
			@timeit to "Initialisaiton" initSimplex = initialisation(prob, SIMPLEX)
			@timeit to "LP Relaxation" UBsimplex, _ = simplex(prob, initSimplex, setvar)
		end 
	end 

	if graphic 
   		# Setup
    	figure("Parametric method and simplex algorithm | "*fname,figsize=(6.5,5))
    	xlabel(L"z^1(x)")
    	ylabel(L"z^2(x)")
    	PyPlot.title("LP Relaxation | "*fname)

		# Show the upper bound set computed by the parametric method
		y_PN11 = [y[1] for y in UBparam.points] 
		y_PN12 = [y[2] for y in UBparam.points]
		scatter(y_PN11, y_PN12, color="green", marker="+", label = "parametric")
		plot(y_PN11, y_PN12, color="green", linewidth=0.75, marker="+",
			markersize=1.0, linestyle=":")

		# Show the upper bound set computed by the simplex algorithm 
		y_PN21 = [y[1] for y in UBsimplex.points] 
		y_PN22 = [y[2] for y in UBsimplex.points]
		scatter(y_PN21, y_PN22, color="red", marker="+", label = "simplex")
		plot(y_PN21, y_PN22, color="red", linewidth=0.75, marker="+",
			markersize=1.0, linestyle=":")

    	legend(bbox_to_anchor=[1,1], loc=0, borderaxespad=0, fontsize = "x-small")
    end 
end

# Calls compareParametric_Simplex on all instances in directory dir 
function testInstancesSimplex(dir::String, graphic=false)

	println("Exemple didactique")
	didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])

	init   = initialisation(didactic, PARAMETRIC_LP)
	setvar = initialSetvar(didactic, init, PARAMETRIC_LP)
	_ = parametricMethod(didactic, init, setvar)

	initSimplex = initialisation(didactic, SIMPLEX)
	setvarSimplex = initialSetvar(didactic, initSimplex, SIMPLEX)
	_ = simplex(didactic, initSimplex, setvarSimplex)
	
	files = readdir(dir)
	for fname in files
		println("\n", fname)
		
		if fname[length(fname)-3:length(fname)] == ".DAT"
			prob = readInstanceMOMKPformatPG(false, dir*fname)
		else
			prob = readInstanceMOMKPformatZL(false, dir*fname)
		end
		compareParametric_Simplex(prob, fname, graphic)
	end	
end

# ----- UBS COMPARISON ------------------------------------------------------- #
# Compare CPU time and memory usage for nIter runs for all methods that compute 
# an upper bound set on the instance in file fname 
function compareUBS(fname::String, nIter::Int)

	didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])

	# Parametric methods 
	init   = initialisation(didactic, PARAMETRIC_LP)
	setvar = initialSetvar(didactic, init, PARAMETRIC_LP)
	_ = parametricMethod(didactic, init, setvar)

	# Dichotomic method 
	initDicho   = initialisation(didactic, DICHOTOMIC)
	setvarDicho = initialSetvar(didactic, initDicho, DICHOTOMIC)
	_ = dichotomicMethod(didactic, initDicho, setvarDicho)
	
	# Simplex algorithm 
	initSimplex   = initialisation(didactic, SIMPLEX) 
	setvarSimplex = initialSetvar(didactic, initSimplex, SIMPLEX)
	_ = simplex(didactic, initSimplex, setvarSimplex)

	if fname[length(fname)-3:length(fname)] == ".DAT"
		prob = readInstanceMOMKPformatPG(false, fname)
	else
		prob = readInstanceMOMKPformatZL(false, fname)
	end

	for iter in 1:nIter 

		# Parametric methods 
		@timeit to "Parametric methods" begin 
			@timeit to "Initialisation" init = initialisation(prob, PARAMETRIC_LP) 
			@timeit to "setvar" setvar = initialSetvar(prob, init, PARAMETRIC_LP)
			@timeit to "LP Relaxation" _ = 
				parametricMethod(prob, init, setvar)
			#=@timeit to "Martello and Toth" _ = 
				martelloAndToth(prob, init, Solution{Float64}(prob))=#
		end 

		# Dichotomic method 
		@timeit to "Dichotomic method" begin 
			@timeit to "Initialisation" initDicho = initialisation(prob, DICHOTOMIC)
			@timeit to "setvar" setvarDicho = initialSetvar(prob, initDicho, DICHOTOMIC)
			@timeit to "Relaxation" _ = 
				dichotomicMethod(prob, initDicho, setvarDicho)
		end 

		# Simplex algorithm 
		@timeit to "Simplex" begin 
			@timeit to "Initialisaiton" initSimplex = initialisation(prob, SIMPLEX)
			@timeit to "setvar" setvarSimplex = initialSetvar(prob, initSimplex, SIMPLEX)
			@timeit to "Relaxation" _ = 
				simplex(prob, initSimplex, setvarSimplex) 
		end 
	end 
end 
 