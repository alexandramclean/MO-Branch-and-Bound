################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Tests : Calcul paramétrique de la relaxation continue                #
################################################################################

include("lpRelaxation.jl")
include("martelloAndToth.jl")
include("dichotomicMethod.jl")
include("simplexAlgorithm.jl")
include("vOptMomkp.jl")
include("parserMomkpPG.jl")
include("parserMomkpZL.jl")
include("displayGraphic.jl")

using TimerOutputs
const to = TimerOutput()

# ----- DICHOTOMIC METHOD ---------------------------------------------------- #
# Produit un graphique affichant la relaxation continue calculée par méthode 
# dichotomique et par méthode paramétrique
function compareLP_DM(name, prob, graphic=false)

	println("Méthode paramétrique")
	println("\tInitialisation : ")
	transpositions, seq, pos = initialisation(prob)
	@timeit to "Relaxation continue" UBparam = parametricMethod(prob, transpositions, seq, pos)
	
	println("Méthode dichotomique")
    UBdicho = dichotomicMethod(prob)

	if graphic 
   		# Setup
    	figure("Test Relaxation Continue | "*name,figsize=(6.5,5))
    	xlabel(L"z^1(x)")
    	ylabel(L"z^2(x)")
    	PyPlot.title("Test Relaxation Continue | "*name)

		# Affichage des points calculés par méthode paramétrique
		y_PN11 = [y[1] for y in UBparam] ; y_PN12 = [y[2] for y in UBparam]
		scatter(y_PN11, y_PN12, color="green", marker="+", label = "parametric")
		plot(y_PN11, y_PN12, color="green", linewidth=0.75, marker="+",
			markersize=1.0, linestyle=":")

		# Affichage des points calculés par méthode dichotomique
		y_PN21 = [y[1] for y in UBdicho] ; y_PN22 = [y[2] for y in UBdicho]
		scatter(y_PN21, y_PN22, color="red", marker="+", label = "dichotomic")
		plot(y_PN21, y_PN22, color="red", linewidth=0.75, marker="+",
			markersize=1.0, linestyle=":")

		#= Affichage des solutions exactes pour le problème non relâché
		if ref != Nothing
			y_N1 = [y[1] for y in ref] ; y_N2 = [y[2] for y in ref]
			scatter(y_N1, y_N2, color="black", marker="+", label = "vOpt")
			plot(y_N1, y_N2, color="black", linewidth=0.75, marker="+",
				markersize=1.0, linestyle=":")
		end=#

    	legend(bbox_to_anchor=[1,1], loc=0, borderaxespad=0, fontsize = "x-small")
    end 
end

# Compare les méthodes paramétrique et dichotomique sur toutes les instances 
# dans le répertoire passé en paramètre 
function testInstances(dir::String, graphic=false) 

	println("Exemple didactique")
	didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
	compareLP_DM("didactic", didactic, graphic)
	
	files = readdir(dir) 
	for fname in files 
		println("\n", fname) 
		if fname[length(fname)-3:length(fname)] == ".DAT"
			prob = readInstanceMOMKPformatPG(false, dir*fname)
		else
			prob = readInstanceMOMKPformatZL(false, dir*fname)
		end
	
		compareLP_DM(fname, prob, graphic)
	end	
end 

# ----- MARTELLO AND TOTH ---------------------------------------------------- #
# Compare CPU times for the LP relaxation and the Martello and Toth upper bound
function compareLP_MT(prob::_MOMKP)

	# Initialisation
	transpositions, seq, pos = initialisation(prob)

	println("LP Relaxation : ")
	@timeit to "LP Relaxation" UBparam = parametricMethod(prob, transpositions, seq, pos)

	# Initialisation
	transpositions, seq, pos = initialisation(prob)
	
	println("Martello and Toth : ")
    @timeit to "Martello and Toth" UB, constraints = martelloAndToth(prob, transpositions, seq, pos)
end 

# Compute LP relaxation and Martello and Toth on all instances in dir
function testInstancesMT(dir::String)

	println("Exemple didactique")
	didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
	compareLP_MT(didactic)
	
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

	@timeit to "Initialisation" transpositions, seq, pos = initialisation(prob)

	var = rand(1:size(prob.P)[2])
	@timeit to "Set variable" newTranspositions, newSeq, newPos = setVariable(transpositions, seq, pos, var)
end

# Calls compareInit_SetVar on all files in directory dir 
function testInstancesSetVar(dir::String)

	println("Exemple didactique")
	didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
	compareInit_SetVar(didactic)
	
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
function compareLP_SP(prob, name, graphic=false)

	#newProb = groupEquivalentItems(prob)

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

	if graphic 
   		# Setup
    	figure("Parametric method and simplex algorithm | "*name,figsize=(6.5,5))
    	xlabel(L"z^1(x)")
    	ylabel(L"z^2(x)")
    	PyPlot.title("LP Relaxation | "*name)

		# Show the upper bound set computed by the parametric method
		y_PN11 = [y[1] for y in UBparam] ; y_PN12 = [y[2] for y in UBparam]
		scatter(y_PN11, y_PN12, color="green", marker="+", label = "parametric")
		plot(y_PN11, y_PN12, color="green", linewidth=0.75, marker="+",
			markersize=1.0, linestyle=":")

		# Show the upper bound set computed by the simplex algorithm 
		y_PN21 = [y[1] for y in UBsimplex] ; y_PN22 = [y[2] for y in UBsimplex]
		scatter(y_PN21, y_PN22, color="red", marker="+", label = "simplex")
		plot(y_PN21, y_PN22, color="red", linewidth=0.75, marker="+",
			markersize=1.0, linestyle=":")

    	legend(bbox_to_anchor=[1,1], loc=0, borderaxespad=0, fontsize = "x-small")
    end 
end

# Calls compareLP_SP on all instances in directory dir 
function testInstancesSimplex(dir::String, graphic=false)

	println("Exemple didactique")
	didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
	transpositions, seq, pos = initialisation(didactic)
	UBparam = parametricMethod(didactic, transpositions, seq, pos)
	seq = simplexInitialisation(didactic)
	UBsimplex = simplex(didactic, seq)
	
	files = readdir(dir)
	for fname in files
		println("\n", fname)
		
		if fname[length(fname)-3:length(fname)] == ".DAT"
			prob = readInstanceMOMKPformatPG(false, dir*fname)
		else
			prob = readInstanceMOMKPformatZL(false, dir*fname)
		end
		compareLP_SP(prob, fname, graphic)
	end
	
end