################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Tests et comparaison de différentes méthodes pour la relaxation      #
################################################################################

include("lpRelaxation.jl")
include("martelloAndToth.jl")
include("dichotomicMethod.jl")
include("simplexAlgorithm.jl")
#include("vOptMomkp.jl")
include("parserMomkpPG.jl")
include("parserMomkpZL.jl")
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
			@timeit to "Initialisation" transpositions, seq, pos = initialisation(prob)
			@timeit to "Relaxation" UBparam = parametricMethod(prob, transpositions, seq, pos)
		end 

		println("Dichotomic method")
		@timeit to "Dichotomic method" begin 
			@timeit to "Initialisation" r1, r2 = utilities(prob)
			@timeit to "Relaxation" UBdicho = dichotomicMethod(prob, r1, r2)
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
		y_PN21 = [y[1] for y in UBdicho] 
		y_PN22 = [y[2] for y in UBdicho]
		scatter(y_PN21, y_PN22, color="red", marker="+", label = "dichotomic")
		plot(y_PN21, y_PN22, color="red", linewidth=0.75, marker="+",
			markersize=1.0, linestyle=":")

		#= Show the non-dominated points
		if ref != Nothing
			y_N1 = [y[1] for y in ref] ; y_N2 = [y[2] for y in ref]
			scatter(y_N1, y_N2, color="black", marker="+", label = "vOpt")
			plot(y_N1, y_N2, color="black", linewidth=0.75, marker="+",
				markersize=1.0, linestyle=":")
		end=#

    	legend(bbox_to_anchor=[1,1], loc=0, borderaxespad=0, fontsize = "x-small")
    end 
end

# Calls compareParametric_Dichotomic on all instances in directory dir 
function testInstances(dir::String, graphic=false) 

	println("Exemple didactique")
	didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
	transpositions, seq, pos = initialisation(didactic)
	_ = parametricMethod(didactic, transpositions, seq, pos)
	r1, r2 = utilities(didactic)
	_ = dichotomicMethod(didactic, r1, r2)
	
	files = readdir(dir) 
	for fname in files 
		println("\n", fname) 
		if fname[length(fname)-3:length(fname)] == ".DAT"
			prob = readInstanceMOMKPformatPG(false, dir*fname)
		else
			prob = readInstanceMOMKPformatZL(false, dir*fname)
		end
	
		compareParametric_Dichotomic(fname, prob, graphic)
	end	
end 

# ----- MARTELLO AND TOTH ---------------------------------------------------- #
# Compare CPU times for the LP relaxation and the Martello and Toth upper bound
function compareLP_MT(prob::_MOMKP)

	@timeit to "\nLP v. MT" begin 
		# Initialisation
		transpositions, seq, pos = initialisation(prob)

		println("LP Relaxation : ")
		@timeit to "LP Relaxation" UBparam = parametricMethod(prob, transpositions, seq, pos)

		# Initialisation
		transpositions, seq, pos = initialisation(prob)
		
		println("Martello and Toth : ")
		@timeit to "Martello and Toth" UB, constraints = martelloAndToth(prob, transpositions, seq, pos)
	end 
end 

# Calls compareLP_MT on all instances in dir
function testInstancesMT(dir::String)

	println("Exemple didactique")
	didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
	transpositions, seq, pos = initialisation(didactic)
	_ = parametricMethod(didactic, transpositions, seq[1:end], pos[1:end]) 
	_ = martelloAndToth(didactic, transpositions, seq[1:end], pos[1:end])
	
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
# Compare CPU times for the initialisation and setVariable functions 
function compareInit_SetVar(prob::_MOMKP)

	@timeit to "\nCompare init and setVar" begin 
		@timeit to "Initialisation" transpositions, seq, pos = initialisation(prob)

		var = rand(1:size(prob.P)[2])
		@timeit to "Set variable" newTranspositions, newSeq, newPos = setVariable(transpositions, seq, pos, var)
	end
end 

# Calls compareInit_SetVar on all files in directory dir 
function testInstancesSetVar(dir::String)

	println("Exemple didactique")
	didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
	transpositions, seq, pos = initialisation(didactic) 
	_ = setVariable(transpositions, seq, pos, rand(1:6))
	
	files = readdir(dir)
	for fname in files
		println("\n", fname)
		
		if fname[length(fname)-3:length(fname)] == ".DAT"
			prob = readInstanceMOMKPformatPG(false, dir*fname)
		else
			prob = readInstanceMOMKPformatZL(false, dir*fname)
		end
		compareInit_SetVar(prob)
	end
	
end

# ----- SIMPLEX ALGORITHM ---------------------------------------------------- #
# Compares the parametric method and the simplex algorithm for computing the 
# LP relaxation of instance prob
# If graphic=true, produces a figure showing both obtained upper bound sets 
function compareParametric_Simplex(prob, name, graphic=false)

	#newProb = groupEquivalentItems(prob)
	@timeit to "\nParametric v. Simplex" begin 
		println("Parametric method")
		@timeit to "Parametric method" begin 
			@timeit to "Initialisation" transpositions, seq, pos = initialisation(prob)
			@timeit to "Relaxation" UBparam = parametricMethod(prob, transpositions, seq, pos)
		end 

		println("Simplex algorithm")
		@timeit to "Simplex algorithm" begin 
			@timeit to "Initialisaiton" seq = simplexInitialisation(prob)
			@timeit to "Relaxation" UBsimplex = simplex(prob, seq)
		end 
	end 

	if graphic 
   		# Setup
    	figure("Parametric method and simplex algorithm | "*name,figsize=(6.5,5))
    	xlabel(L"z^1(x)")
    	ylabel(L"z^2(x)")
    	PyPlot.title("LP Relaxation | "*name)

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
	transpositions, seq, pos = initialisation(didactic)
	_ = parametricMethod(didactic, transpositions, seq, pos)
	seq = simplexInitialisation(didactic)
	_ = simplex(didactic, seq)
	
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
function compareUBS(fname::String, nIter::Int)

	didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])

	# Parametric methods 
	transpositions, seq, pos = initialisation(didactic)
	_ = parametricMethod(didactic, transpositions, seq[1:end], pos[1:end])
	_ = martelloAndToth(didactic, transpositions, seq[1:end], pos[1:end]) 

	# Dichotomic method 
	r1, r2 = utilities(didactic)
	_ = dichotomicMethod(didactic, r1, r2)
	
	# Simplex algorithm 
	seq = simplexInitialisation(didactic) 
	_ = simplex(didactic, seq)

	if fname[length(fname)-3:length(fname)] == ".DAT"
		prob = readInstanceMOMKPformatPG(false, fname)
	else
		prob = readInstanceMOMKPformatZL(false, fname)
	end

	for iter in 1:nIter 

		# Parametric methods 
		@timeit to "Parametric methods" begin 
			@timeit to "Initialisation" transpositions, seq, pos = initialisation(prob) 
			@timeit to "LP Relaxation" _ = parametricMethod(prob, transpositions, seq[1:end], pos[1:end])
			@timeit to "Martello and Toth" _ = martelloAndToth(prob, transpositions, seq[1:end], pos[1:end])
		end 

		# Dichotomic method 
		@timeit to "Dichotomic method" begin 
			@timeit to "Initialisation" r1, r2 = utilities(prob) 
			@timeit to "Relaxation" _ = dichotomicMethod(prob, r1, r2)
		end 

		# Simplex algorithm 
		@timeit to "Simplex" begin 
			@timeit to "Initialisaiton" seq = simplexInitialisation(prob)
			@timeit to "Relaxation" _ = simplex(prob, seq) 
		end 
	end 
end 