################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Tests : Calcul paramétrique de la relaxation continue                #
################################################################################

include("relaxationContinue.jl")
include("dichotomicMethod.jl")
include("vOptMomkp.jl")
include("parserMomkpPG.jl")
include("parserMomkpZL.jl")
include("displayGraphic.jl")

# Produit un graphique affichant l'ensemble des points non-dominés pour un 
# problème donné, ainsi que la relaxation calculée par méthode dichotomique et 
# par méthode paramétrique 
function testComparaison(name, prob, ref)
	
	@time UBparam = relaxationContinue(prob)
    @time UBdicho = dichotomicMethod(prob) 

    # Setup
    figure("Test",figsize=(6.5,5))
    xlabel(L"z^1(x)")
    ylabel(L"z^2(x)")
    PyPlot.title("Test Relaxation Continue | "*name)

	# Affichage des points calculés par méthode paramétrique
	y_PN11 = [y.z[1] for y in UBparam] ; y_PN12 = [y.z[2] for y in UBparam]
    scatter(y_PN11, y_PN12, color="green", marker="+", label = "parametric")
    plot(y_PN11, y_PN12, color="green", linewidth=0.75, marker="+", 
    	markersize=1.0, linestyle=":")
    
    # Affichage des points calculés par méthode dichotomique
    y_PN21 = [y[1] for y in UBdicho] ; y_PN22 = [y[2] for y in UBdicho]
    scatter(y_PN21, y_PN22, color="red", marker="+", label = "dichotomic")
    plot(y_PN21, y_PN22, color="red", linewidth=0.75, marker="+", 
    	markersize=1.0, linestyle=":")

	# Affichage des solutions exactes pour le problème non relâché
	y_N1 = [y[1] for y in ref] ; y_N2 = [y[2] for y in ref]
    scatter(y_N1, y_N2, color="black", marker="+", label = "vOpt")
    plot(y_N1, y_N2, color="black", linewidth=0.75, marker="+", 
    	markersize=1.0, linestyle=":")

    legend(bbox_to_anchor=[1,1], loc=0, borderaxespad=0, fontsize = "x-small")
end

#= Exemple didactique
prob = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
ref, _ = vSolveBi01IP(GLPK.Optimizer, didactic.P, didactic.W, didactic.ω)
testComparaison("didactic", didactic, ref) =#

# ----------
fname = "../instancesPG/set1/ZL28.DAT"
ref = [[902.0, 1193.0], [962.0, 1187.0], [967.0, 1186.0], [994.0, 1163.0], [1080.0, 1155.0], [1108.0, 1131.0], [1117.0, 1122.0], [1119.0, 1118.0], [1142.0, 1083.0], [1148.0, 1056.0], [1159.0, 1044.0], [1167.0, 1005.0], [1180.0, 1004.0], [1181.0, 1003.0], [1182.0, 979.0], [1198.0, 977.0], [1219.0, 969.0], [1220.0, 818.0]]

if fname[length(fname)-3:length(fname)] == ".DAT"
    prob = readInstanceMOMKPformatPG(false, fname)
else
    prob = readInstanceMOMKPformatZL(false, fname)
end

testComparaison(fname, prob, ref) 






