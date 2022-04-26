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

# ----- REFERENCE SETS ------------------------------------------------------- #
# Compute reference sets for all files in the specified folder 
function computeReferenceSets(dir::String)

	files = readdir(dir*"dat/")

	for fname in files 

		println("\n", fname)

		# Read the instance in the file 
		if fname[length(fname)-3:length(fname)] == ".DAT"
			prob = readInstanceMOMKPformatPG(false, dir*"dat/"*fname)
		else
			prob = readInstanceMOMKPformatZL(false, dir*"dat/"*fname)
		end

		# Transform the multi-dimensional problem into a mono-dimensional problem
		prob = multiToMonoDimensional(prob)

		# Solve the problem with vOpt
		start = time() 
		ref, _ = vSolveBi01IP(GLPK.Optimizer, prob.P, prob.W, prob.ω)
		elapsed = time() - start 

		# Write in a file 
		open(dir*"res/ref_"*fname, "w") do io 
			write(io, fname*"\n")
			write(io, string(elapsed))
			for y in ref 
				write(io, "\n"*string(y[1])*" "*string(y[2]))
			end 
		end; 
	end 
end 

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
			Ldicho = Vector{Solution{Rational{Int}}}()
			@timeit to "Relaxation" UBdicho = 
				dichotomicMethod(prob, Ldicho, initDicho)
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
		y_PN21 = [y.point[1] for y in UBdicho] 
		y_PN22 = [y.point[2] for y in UBdicho]
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
	didactic  = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])

	L = PrimalBoundSet{Float64}()
	init = initialisation(didactic, PARAMETRIC_LP)
	setvar = initialSetvar(didactic, init, PARAMETRIC_LP)
	_ = parametricMethod(didactic, init, setvar)

	L = Vector{Solution{Rational{Int}}}()
	initDicho = initialisation(didactic, DICHOTOMIC)
	_ = dichotomicMethod(didactic, L, initDicho)
	
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
# Compare CPU times for the initialisation and setVariable functions 
function compareInit_SetVar(prob::_MOMKP)

	@timeit to "\nCompare init and setVar" begin 
		@timeit to "Initialisation" init = initialisation(prob, PARAMETRIC_LP)

		var = rand(1:size(prob.P)[2])
		@timeit to "Set variable" newInit = setVariable(init, var, PARAMETRIC_LP)
	end
end 

# Calls compareInit_SetVar on all files in directory dir 
function testInstancesInit_SetVar(dir::String)

	println("Exemple didactique")
	didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
	init = initialisation(didactic, PARAMETRIC_LP) 
	_ = setVariable(init, rand(1:6), PARAMETRIC_LP)
	
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

# Compute the upper bound set with a variable set to 0 or 1 with the parametric,
# dichotomic, and simplex methods 
function compareSetVar(fname::String, 
					   var::Union{Int,Nothing} = nothing, 
					   val::Union{Int,Nothing} = nothing)

	println(basename(fname))

	# Read the instance in the file  
	if fname[length(fname)-3:length(fname)] == ".DAT"
		prob = readInstanceMOMKPformatPG(false, fname)
	else
		prob = readInstanceMOMKPformatZL(false, fname)
	end

	n = size(prob.P)[2]

	if var === nothing
		# Generate a random variable 
		var = rand(1:n)
	end 

	if val === nothing 
		# Generate a random value 
		val = rand(0:1)
	end 

	@assert var >= 1 && var <= n "The variable you are trying to set doesn't exist"
	@assert val == 0 || val == 1 "Variables can only be set to 0 or 1"

	println("(", var, ",", val, ")")
	# Initialisaiton
	init = initialisation(prob, PARAMETRIC_LP)

	# Variable setting 
	newInit = setVariable(init, var, PARAMETRIC_LP)

	# Parametric method 
	L = PrimalBoundSet{Float64}() 

	solInit = Solution{Float64}(prob)
	solInit.X[var] = val 
	solInit.z     += val*prob.P[:,var]
	if val == 1 
		solInit.ω_ -= prob.W[1,var]
	end 

	UBparam = parametricMethod(prob, L, newInit, solInit)

	# Dichotomic method 
	L = PrimalBoundSet{Rational{Int}}() 

	# Variable setting 
	newInit = setVariable(init, var, DICHOTOMIC)	

	solInit = Solution{Rational{Int}}(prob)
	solInit.X[var] = val 
	solInit.z     += val*prob.P[:,var]
	if val == 1 
		solInit.ω_ -= prob.W[1,var]
	end 

	UBdicho = dichotomicMethod(prob, L, newInit, solInit)	

	# Simplex method 
	L = PrimalBoundSet{Float64}()

	# Variable setting 
	newInit = setVariable(init, var, SIMPLEX)

	solInit = Solution{Float64}(prob)
	solInit.X[var] = val 
	solInit.z     += val*prob.P[:,var]
	if val == 1 
		solInit.ω_ -= prob.W[1,var]
	end 

	UBsimplex = simplex(prob, L, newInit, solInit)
	
	# Setup
	figure("Variable setting | "*basename(fname),figsize=(6.5,5))
	xlabel(L"z^1(x)")
	ylabel(L"z^2(x)")
	PyPlot.title("Parametric, dichotomic, and simplex methods | "*basename(fname))

	# Show the points computed by the parametric method 
	y_PN11 = [y[1] for y in UBparam.points] 
	y_PN12 = [y[2] for y in UBparam.points]
	scatter(y_PN11, y_PN12, color="black", marker="o", label = "parametric")
	plot(y_PN11, y_PN12, color="black", linewidth=0.75, marker="o",
		markersize=1.0, linestyle=":")

	# Show the points computed by the dichotomic method 
	y_PN21 = [y[1] for y in UBdicho.points] 
	y_PN22 = [y[2] for y in UBdicho.points]
	scatter(y_PN21, y_PN22, color="red", marker="*", label = "dichotomic")
	plot(y_PN21, y_PN22, color="red", linewidth=0.75, marker="*",
		markersize=1.0, linestyle=":")

	# Show the points computed by the simplex algorithm 
	y_PN31 = [y[1] for y in UBsimplex.points] 
	y_PN32 = [y[2] for y in UBsimplex.points]
	scatter(y_PN31, y_PN32, color="green", marker="+", label = "simplex")
	plot(y_PN31, y_PN32, color="green", linewidth=0.75, marker="+",
		markersize=1.0, linestyle=":")


	legend(bbox_to_anchor=[1,1], loc=0, borderaxespad=0, fontsize = "x-small")
end 

function testInstancesSetVar(dir::String)

	files = readdir(dir)

	for fname in files 
		compareSetVar(dir*fname)
	end 
end 

# ----- SIMPLEX ALGORITHM ---------------------------------------------------- #
# Compares the parametric method and the simplex algorithm for computing the 
# LP relaxation of instance prob
# If graphic=true, produces a figure showing both obtained upper bound sets 
function compareParametric_Simplex(prob::_MOMKP, fname::String, graphic=false)

	#newProb = groupEquivalentItems(prob)
	@timeit to "\nParametric v. Simplex" begin 
		println("Parametric method")
		@timeit to "Parametric method" begin 
			@timeit to "Initialisation" init = initialisation(prob, PARAMETRIC_LP)
			L = PrimalBoundSet{Float64}()
			@timeit to "Relaxation" UBparam = 
				parametricMethod(prob, L, init, Solution{Float64}(prob))
		end 

		println("Simplex algorithm")
		@timeit to "Simplex algorithm" begin 
			@timeit to "Initialisaiton" initSimplex = initialisation(prob, SIMPLEX)
			L = PrimalBoundSet{Float64}()
			@timeit to "Relaxation" UBsimplex = 
				simplex(prob, L, initSimplex, Solution{Float64}(prob))
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

	L = PrimalBoundSet{Float64}()
	init = initialisation(didactic, PARAMETRIC_LP)
	_ = parametricMethod(didactic, L, init, Solution{Float64}(didactic))

	L = PrimalBoundSet{Float64}()
	initSimplex = initialisation(didactic, SIMPLEX)
	_ = simplex(didactic, L, initSimplex, Solution{Float64}(didactic))
	
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
	init = initialisation(didactic, PARAMETRIC_LP)
	L = PrimalBoundSet{Float64}()
	_ = parametricMethod(didactic, L, init, Solution{Float64}(didactic))
	_ = martelloAndToth(didactic, init, Solution{Float64}(didactic)) 

	# Dichotomic method 
	initDicho = initialisation(didactic, DICHOTOMIC)
	L = PrimalBoundSet{Rational{Int}}()
	_ = dichotomicMethod(didactic, L, initDicho, Solution{Rational{Int}}(didactic))
	
	# Simplex algorithm 
	initSimplex = initialisation(didactic, SIMPLEX) 
	L = PrimalBoundSet{Float64}()
	_ = simplex(didactic, L, initSimplex, Solution{Float64}(didactic))

	if fname[length(fname)-3:length(fname)] == ".DAT"
		prob = readInstanceMOMKPformatPG(false, fname)
	else
		prob = readInstanceMOMKPformatZL(false, fname)
	end

	for iter in 1:nIter 

		# Parametric methods 
		L = PrimalBoundSet{Float64}() 
		@timeit to "Parametric methods" begin 
			@timeit to "Initialisation" init = initialisation(prob, PARAMETRIC_LP) 
			@timeit to "setvar" setvar = initialSetvar(prob, init, PARAMETRIC_LP)
			@timeit to "LP Relaxation" _ = 
				parametricMethod(prob, L, init, setvar, false)
			#=@timeit to "Martello and Toth" _ = 
				martelloAndToth(prob, init, Solution{Float64}(prob))=#
		end 

		# Dichotomic method 
		L = PrimalBoundSet{Rational{Int}}()
		@timeit to "Dichotomic method" begin 
			@timeit to "Initialisation" initDicho = initialisation(prob, DICHOTOMIC)
			@timeit to "Relaxation" _ = 
				dichotomicMethod(prob, L, initDicho, Solution{Rational{Int}}(prob))
		end 

		# Simplex algorithm 
		L = PrimalBoundSet{Float64}()
		@timeit to "Simplex" begin 
			@timeit to "Initialisaiton" initSimplex = initialisation(prob, SIMPLEX)
			@timeit to "Relaxation" _ = 
				simplex(prob, L, initSimplex, Solution{Float64}(prob)) 
		end 
	end 
end 

# ----- BRANCH-AND-BOUND ----------------------------------------------------- #
# Tests the branch-and-bound algorithm on the instance in file fname 
function testBranchAndBound(fname::String, 
							ref::Vector{Vector{Float64}},
							method::Method=PARAMETRIC_LP,
							interrupt::Bool=false,
							initialisation::Bool=true) where T<:Real 

	println(basename(fname))

	# Read the instance in the file 
	if fname[length(fname)-3:length(fname)] == ".DAT"
		prob = readInstanceMOMKPformatPG(false, fname)
	else
		prob = readInstanceMOMKPformatZL(false, fname)
	end

	if initialisation 
		#@timeit to "Supported efficient" L = dichotomicMethod(prob) 
		L = PrimalBoundSet(ref) 
	else 
		if method == DICHOTOMIC 
			L = PrimalBoundSet{Rational{Int}}()
		else 
			L = PrimalBoundSet{Float64}()
		end 
	end 

	# Branch-and-bound 
	@timeit to "Branch-and-bound" L = 
		branchAndBound(prob, L, method, interrupt)

	# Plot the obtained set of solutions 
	plotYN(basename(fname), ref, L.solutions)

	@assert ref == [sol.z for sol in L.solutions] "The solutions obtained by the 
	branch-and-bound algorithm must be identical to those in the reference set"
end 

# Tests the branch-and-bound algorithm on all instances in directory dir 
function testInstancesBranchAndBound(dir::String, 
									 method::Method=PARAMETRIC_LP,
									 interrupt::Bool=false,
									 initialisation::Bool=true)

	# Exemple didactique 
	didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
	if initialisation 
		L = dichotomicMethod(didactic) 
	else 
		if method == DICHOTOMIC 
			L = PrimalBoundSet{Rational{Int}}()
		else 
			L = PrimalBoundSet{Float64}()
		end 
	end 

	# Branch-and-bound 
	L = branchAndBound(didactic, L, method, interrupt)

	
	files = readdir(dir*"dat/")
	for fname in files 
		# Get the reference set 
		ref = readReferenceSet(dir*"res/ref_"*fname)

		testBranchAndBound(dir*"dat/"*fname, ref, method, interrupt, initialisation)
	end 

end 