################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Test functions for the branch-and-bound algorithm                    #
################################################################################

include("lpRelaxation.jl")
include("martelloAndToth.jl")
include("dichotomicMethod.jl")
include("simplexAlgorithm.jl")
include("branch-and-bound.jl")
include("efficientSupportedSolutions.jl")
include("vOptMomkp.jl")
include("parsers.jl")

using TimerOutputs
const to = TimerOutput()

# ----- REFERENCE SETS ------------------------------------------------------- #
# Compute reference sets for all files in the specified folder and stores them 
# in a format readable by the readReferenceSet function 
function computeReferenceSets(dir::String, groupEquivItems::Bool=false)

	files = readdir(dir*"dat/")

	for fname in files 

		println("\n", fname)

		# Read the instance in the file 
		prob = readInstance(dir*"dat/"*fname)

		# Transform the multi-dimensional problem into a mono-dimensional problem
		prob = multiToMonoDimensional(prob)
		groupEquivItems ? prob = groupEquivalentItems(prob) : nothing 

		# Solve the problem with vOpt
		start = time() 
		ref, _ = vSolveBi01IP(GLPK.Optimizer, prob.P, prob.W, prob.ω)
		elapsed = time() - start 

		# Write in a file 
		if groupEquivItems
			resFile = dir*"resGrouped/ref_"*fname 
		else 
			resFile = dir*"res/ref_"*fname 
		end 

		open(groupEquivItems, "w") do io 
			write(io, fname*"\n")
			write(io, string(elapsed))
			for y in ref 
				write(io, "\n"*string(y[1])*" "*string(y[2]))
			end 
		end; 
	end 
end 

# ----- BRANCH-AND-BOUND ----------------------------------------------------- #
# Tests the branch-and-bound algorithm on the instance in file fname 
function testBranchAndBound(fname::String, 
							ref::Vector{Vector{Float64}},
							method::Method=PARAMETRIC_LP,
							initialisation::Bool=false,
							interrupt::Bool=false							
						   ) where T<:Real 

	println("\n", basename(fname))

	# Read the instance in the file 
	prob = readInstance(fname)

	#groupEquivItems ? prob = groupEquivalentItems(prob) : nothing 

	if initialisation 
		@timeit to "Supported efficient" L = dichotomicMethod(prob) 
	else 
		if method == DICHOTOMIC 
			L = Vector{Solution{Rational{Int}}}()
		else 
			L = Vector{Solution{Float64}}()
		end 
	end 

	# Branch-and-bound 
	@timeit to "Branch-and-bound" L = branchAndBound(prob, L, method, interrupt)

	if !(ref == [sol.z for sol in L])
		println("The solutions obtained by the branch-and-bound algorithm must be identical to those in the reference set")
	end
end 

# Tests the branch-and-bound algorithm on all instances in directory dir 
function testInstancesBranchAndBound(dir::String, 
									 method::Method=PARAMETRIC_LP,
									 initialisation::Bool=false,
									 interrupt::Bool=false)

	# Exemple didactique 
	didactic = didacticInstance()
	if initialisation 
		L = dichotomicMethod(didactic) 
	else 
		if method == DICHOTOMIC 
			L = Vector{Solution{Rational{Int}}}()
		else 
			L = Vector{Solution{Float64}}()
		end 
	end 

	# Branch-and-bound 
	L = branchAndBound(didactic, L, method, interrupt)

	files = readdir(dir*"dat/")
	for fname in files 
		# Get the reference set 
		#if groupEquivItems
		#	ref = readReferenceSet(dir*"resGrouped/ref_"*fname)
		#else 
			ref = readReferenceSet(dir*"res/ref_"*fname)
		#end 

		testBranchAndBound(dir*"dat/"*fname, ref, method, initialisation, interrupt)
	end 
end 

# Test the branch-and-bound algorithm on one instance for nIter runs 
function testBranchAndBound(fname::String, # Instance 
							# Method used to compute upper bound sets 
							method::Method=PARAMETRIC_LP,
							# Supported efficient solutions
							initialisation::Bool=false,
							# Whether to interrupt the computation of the UBS 
							interrupt::Bool=false,
							nIter::Int=10)

	# Small test problem  
	didactic = didacticInstance()
	if initialisation 
		L = dichotomicMethod(didactic) 
	else 
		if method == DICHOTOMIC 
			L = Vector{Solution{Rational{Int}}}()
		else 
			L = Vector{Solution{Float64}}()
		end 
	end 

	# Branch-and-bound 
	L = branchAndBound(didactic, L, method, interrupt)

	prob = readInstance(fname)

	for iter in 1:nIter 

		if initialisation 
			@timeit to "Supported efficient" L = 
				dichotomicMethod(prob) 
		else 
			if method == DICHOTOMIC 
				L = Vector{Solution{Rational{Int}}}()
			else 
				L = Vector{Solution{Float64}}()
			end 
		end

		# Branch-and-bound 
		@timeit to "Branch-and-bound" L = 
			branchAndBound(prob, L, method, interrupt)
	end 
end 