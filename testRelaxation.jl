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
didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
ref, _ = vSolveBi01IP(GLPK.Optimizer, didactic.P, didactic.W, didactic.ω)
testComparaison("didactic", didactic, ref) =#

# ----------
fname = "ZL105.DAT"
ref = [[3190.0, 4307.0], [3191.0, 4302.0], [3231.0, 4301.0], [3241.0, 4300.0], [3335.0, 4298.0], [3339.0, 4295.0], [3353.0, 4294.0], [3355.0, 4292.0], [3405.0, 4291.0], [3407.0, 4289.0], [3410.0, 4284.0], [3414.0, 4283.0], [3419.0, 4282.0], [3453.0, 4281.0], [3483.0, 4279.0], [3493.0, 4278.0], [3546.0, 4274.0], [3558.0, 4267.0], [3570.0, 4261.0], [3588.0, 4260.0], [3591.0, 4253.0], [3626.0, 4252.0], [3640.0, 4248.0], [3675.0, 4241.0], [3684.0, 4231.0], [3701.0, 4225.0], [3716.0, 4220.0], [3717.0, 4219.0], [3724.0, 4214.0], [3745.0, 4212.0], [3751.0, 4208.0], [3776.0, 4206.0], [3786.0, 4193.0], [3787.0, 4191.0], [3807.0, 4189.0], [3811.0, 4177.0], [3817.0, 4176.0], [3825.0, 4170.0], [3826.0, 4169.0], [3835.0, 4167.0], [3860.0, 4163.0], [3862.0, 4152.0], [3873.0, 4149.0], [3887.0, 4144.0], [3893.0, 4135.0], [3897.0, 4131.0], [3918.0, 4127.0], [3928.0, 4114.0], [3936.0, 4108.0], [3946.0, 4095.0], [3948.0, 4088.0], [3955.0, 4082.0], [3964.0, 4080.0], [3966.0, 4079.0], [3968.0, 4072.0], [3969.0, 4069.0], [3974.0, 4067.0], [3977.0, 4064.0], [3978.0, 4062.0], [3984.0, 4060.0], [3986.0, 4057.0], [3999.0, 4055.0], [4009.0, 4045.0], [4014.0, 4035.0], [4024.0, 4024.0], [4025.0, 4013.0], [4034.0, 4011.0], [4035.0, 4002.0], [4038.0, 4001.0], [4046.0, 4000.0], [4047.0, 3993.0], [4048.0, 3990.0], [4053.0, 3982.0], [4057.0, 3978.0], [4062.0, 3971.0], [4067.0, 3965.0], [4069.0, 3959.0], [4078.0, 3958.0], [4082.0, 3954.0], [4086.0, 3947.0], [4091.0, 3943.0], [4092.0, 3933.0], [4116.0, 3928.0], [4119.0, 3909.0], [4124.0, 3905.0], [4128.0, 3898.0], [4133.0, 3891.0], [4139.0, 3886.0], [4140.0, 3873.0], [4145.0, 3867.0], [4150.0, 3860.0], [4159.0, 3850.0], [4165.0, 3844.0], [4168.0, 3839.0], [4171.0, 3825.0], [4174.0, 3824.0], [4175.0, 3822.0], [4179.0, 3814.0], [4185.0, 3808.0], [4190.0, 3802.0], [4191.0, 3793.0], [4199.0, 3792.0], [4208.0, 3781.0], [4215.0, 3763.0], [4216.0, 3760.0], [4227.0, 3749.0], [4231.0, 3728.0], [4232.0, 3715.0], [4233.0, 3709.0], [4236.0, 3700.0], [4241.0, 3698.0], [4249.0, 3696.0], [4252.0, 3679.0], [4253.0, 3672.0], [4256.0, 3664.0], [4258.0, 3660.0], [4261.0, 3654.0], [4264.0, 3639.0], [4266.0, 3629.0], [4274.0, 3625.0], [4276.0, 3612.0], [4277.0, 3598.0], [4281.0, 3596.0], [4283.0, 3593.0], [4285.0, 3573.0], [4290.0, 3563.0], [4293.0, 3558.0], [4294.0, 3548.0], [4298.0, 3535.0], [4301.0, 3520.0], [4302.0, 3517.0], [4303.0, 3507.0], [4306.0, 3502.0], [4308.0, 3499.0], [4311.0, 3497.0], [4320.0, 3487.0], [4323.0, 3467.0], [4326.0, 3434.0], [4327.0, 3426.0], [4330.0, 3422.0], [4331.0, 3397.0], [4336.0, 3390.0], [4337.0, 3374.0], [4338.0, 3364.0], [4341.0, 3362.0], [4344.0, 3325.0], [4352.0, 3321.0], [4354.0, 3302.0], [4357.0, 3282.0], [4358.0, 3277.0], [4359.0, 3203.0], [4368.0, 3193.0]]


if fname[length(fname)-3:length(fname)] == ".DAT"
    prob = readInstanceMOMKPformatPG(false, fname)
else
    prob = readInstanceMOMKPformatZL(false, fname)
end

testComparaison(fname, prob, ref)
