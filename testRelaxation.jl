################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Tests : Calcul paramétrique de la relaxation continue                #
################################################################################

include("lpRelaxation.jl")
include("martelloAndToth.jl")
include("dichotomicMethod.jl")
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

# Exemple didactique
function testDidactic(graphic=false)
	didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
	compareLP_DM("didactic", didactic, graphic)
end

# Test sur une instance contenue dans le fichier fname
function testFile(fname::String, graphic=false)

	if fname[length(fname)-3:length(fname)] == ".DAT"
    	prob = readInstanceMOMKPformatPG(false, fname)
	else
    	prob = readInstanceMOMKPformatZL(false, fname)
	end

	compareLP_DM(fname, prob, graphic)
end

# Compare les méthodes paramétrique et dichotomique sur toutes les instances 
function testInstances(dir::String, graphic=false) 

	println("Exemple didactique")
	testDidactic(graphic)
	
	files = readdir(dir) 
	for fname in files 
		println("\n", fname) 
		testFile(dir*fname, graphic) 
	end	
end 

# ----- MARTELLO AND TOTH ---------------------------------------------------- #
# Compare execution times for the LP relaxation and Martello and Toth
function compareLP_MT(prob::_MOMKP)

	# Initialisation
	transpositions, seq, pos = initialisation(prob)

	println("Relaxation continue : ")
	@time UBparam = parametricMethod(prob, transpositions, seq, pos)

	# Initialisation
	transpositions, seq, pos = initialisation(prob)
	
	println("Martello et Toth : ")
    @time UB, constraints = martelloAndToth(prob, transpositions, seq, pos)
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
# Compare execution times for the initialisation and setVariable functions 
function compareLP_SetVar(prob::_MOMKP)

	# Initialisation
	#println("Initialisation : ")
	transpositions, seq, pos = initialisation(prob)

	#println("Set variable : ")
	var = 1
	@timeit to "Set variable" newTranspositions, newSeq, newPos = setVariable(transpositions, seq, pos, var)
	
end

function testInstancesSetVar(dir::String)

	println("Exemple didactique")
	didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
	compareLP_SetVar(didactic)
	
	files = readdir(dir)
	for fname in files
		println("\n", fname)
		
		if fname[length(fname)-3:length(fname)] == ".DAT"
			prob = readInstanceMOMKPformatPG(false, dir*fname)
		else
			prob = readInstanceMOMKPformatZL(false, dir*fname)
		end
		compareLP_SetVar(prob)
	end
	
end

# ----- SIMPLEX ALGORITHM ---------------------------------------------------- #
function compareLP_SP(name, prob, graphic=false)

	println("Méthode paramétrique")
	println("\tInitialisation : ")
	@time transpositions, seq, pos = initialisation(prob)
	@time UBparam = parametricMethod(prob, transpositions, seq, pos)
	
	println("Algorithme du simplexe")
    @time UBsimplexe = dichotomicMethod(prob)

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
		y_PN21 = [y[1] for y in UBsimplexe] ; y_PN22 = [y[2] for y in UBsimplexe]
		scatter(y_PN21, y_PN22, color="red", marker="+", label = "simplex")
		plot(y_PN21, y_PN22, color="red", linewidth=0.75, marker="+",
			markersize=1.0, linestyle=":")

    	legend(bbox_to_anchor=[1,1], loc=0, borderaxespad=0, fontsize = "x-small")
    end 
end