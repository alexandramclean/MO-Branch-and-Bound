################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Tests : Calcul paramétrique de la relaxation continue                #
################################################################################

include("relaxationContinue.jl")
include("vOptMomkp.jl")
include("parserMomkpPG.jl")
include("parserMomkpZL.jl")
include("displayGraphic.jl")

# Graphique avec l'ensemble bornant supérieur et l'ensemble des points non-dominés
function testRelaxation(name, prob, ref)
	
    upperBound = relaxationContinue(prob) 
	#=
    for sol in upperBound 
    	println(sol.z) 
    end =#

    # Setup
    figure("Test",figsize=(6.5,5))
    xlabel(L"z^1(x)")
    ylabel(L"z^2(x)")
    PyPlot.title("Test Relaxation Continue | "*name)

    y_PN1 = [x.z[1] for x in upperBound]
    y_PN2 = [x.z[2] for x in upperBound]
    y_N1 = [x[1] for x in ref]
    y_N2 = [x[2] for x in ref]
    
    # display only the points corresponding to non-dominated points
    scatter(y_PN1, y_PN2, color="green", marker="+", label = "UBS")
    # display segments joining adjacent non-dominated points
    plot(y_PN1, y_PN2, color="green", linewidth=0.75, marker="+", markersize=1.0, linestyle=":")

    scatter(y_N1, y_N2, color="black", marker="+", label = "vOpt")
    # display segments joining non-dominated points and their corners points
    Env1,Env2 = computeCornerPointsLowerEnvelop(y_N1, y_N2)
    plot(Env1, Env2, color="black", linewidth=0.75, marker="+", markersize=1.0, linestyle=":")

    legend(bbox_to_anchor=[1,1], loc=0, borderaxespad=0, fontsize = "x-small")
end

#= Exemple didactique
didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
ref, _ = vSolveBi01IP(GLPK.Optimizer, didactic.P, didactic.W, didactic.ω)
@time testRelaxation("didactic", didactic, ref) =#

# ----------
fname = "../instancesPG/set1/ZL100.DAT"
#ref = [[902.0, 1193.0], [962.0, 1187.0], [967.0, 1186.0], [994.0, 1163.0], [1080.0, 1155.0], [1108.0, 1131.0], [1117.0, 1122.0], [1119.0, 1118.0], [1142.0, 1083.0], [1148.0, 1056.0], [1159.0, 1044.0], [1167.0, 1005.0], [1180.0, 1004.0], [1181.0, 1003.0], [1182.0, 979.0], [1198.0, 977.0], [1219.0, 969.0], [1220.0, 818.0]]

ref = [[3235.0, 4037.0], [3246.0, 4031.0], [3271.0, 4029.0], [3340.0, 4028.0], [3350.0, 4024.0], [3364.0, 4021.0], [3415.0, 4017.0], [3431.0, 4009.0], [3453.0, 4003.0], [3466.0, 4001.0], [3469.0, 4000.0], [3481.0, 3998.0], [3485.0, 3993.0], [3490.0, 3990.0], [3497.0, 3989.0], [3511.0, 3987.0], [3512.0, 3986.0], [3540.0, 3982.0], [3556.0, 3975.0], [3568.0, 3971.0], [3582.0, 3970.0], [3595.0, 3966.0], [3616.0, 3956.0], [3632.0, 3949.0], [3635.0, 3938.0], [3645.0, 3937.0], [3653.0, 3936.0], [3660.0, 3935.0], [3666.0, 3926.0], [3672.0, 3925.0], [3688.0, 3919.0], [3690.0, 3914.0], [3704.0, 3912.0], [3732.0, 3908.0], [3749.0, 3893.0], [3751.0, 3891.0], [3753.0, 3885.0], [3762.0, 3881.0], [3763.0, 3878.0], [3778.0, 3875.0], [3779.0, 3871.0], [3787.0, 3870.0], [3791.0, 3863.0], [3797.0, 3862.0], [3806.0, 3860.0], [3822.0, 3857.0], [3825.0, 3851.0], [3831.0, 3848.0], [3832.0, 3845.0], [3833.0, 3843.0], [3834.0, 3842.0], [3838.0, 3839.0], [3839.0, 3836.0], [3855.0, 3832.0], [3860.0, 3830.0], [3870.0, 3826.0], [3879.0, 3817.0], [3884.0, 3813.0], [3888.0, 3812.0], [3907.0, 3806.0], [3909.0, 3801.0], [3918.0, 3792.0], [3924.0, 3786.0], [3927.0, 3783.0], [3929.0, 3778.0], [3938.0, 3773.0], [3952.0, 3769.0], [3963.0, 3761.0], [3964.0, 3748.0], [3969.0, 3747.0], [3970.0, 3739.0], [3977.0, 3738.0], [3990.0, 3731.0], [3994.0, 3724.0], [3996.0, 3722.0], [3997.0, 3718.0], [4009.0, 3717.0], [4010.0, 3716.0], [4014.0, 3705.0], [4019.0, 3700.0], [4041.0, 3697.0], [4050.0, 3680.0], [4056.0, 3666.0], [4064.0, 3660.0], [4071.0, 3649.0], [4074.0, 3646.0], [4082.0, 3640.0], [4089.0, 3635.0], [4094.0, 3612.0], [4098.0, 3605.0], [4100.0, 3603.0], [4102.0, 3602.0], [4121.0, 3595.0], [4122.0, 3569.0], [4125.0, 3566.0], [4128.0, 3560.0], [4130.0, 3555.0], [4136.0, 3554.0], [4139.0, 3546.0], [4142.0, 3537.0], [4151.0, 3534.0], [4152.0, 3528.0], [4164.0, 3521.0], [4168.0, 3505.0], [4174.0, 3499.0], [4175.0, 3491.0], [4182.0, 3478.0], [4185.0, 3466.0], [4197.0, 3462.0], [4199.0, 3443.0], [4205.0, 3437.0], [4215.0, 3423.0], [4218.0, 3400.0], [4220.0, 3368.0], [4230.0, 3367.0], [4236.0, 3338.0], [4246.0, 3319.0], [4248.0, 3300.0], [4250.0, 3278.0], [4262.0, 3274.0], [4266.0, 3215.0]]

if fname[length(fname)-3:length(fname)] == ".DAT"
    prob = readInstanceMOMKPformatPG(false, fname)
else
    prob = readInstanceMOMKPformatZL(false, fname)
end

@time testRelaxation(fname, prob, ref) 
