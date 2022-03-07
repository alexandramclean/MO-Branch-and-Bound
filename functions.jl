################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Fonctions auxiliaires                                                #
################################################################################

include("dataStructures.jl")

# ----- DOMINANCE ------------------------------------------------------------ #
# Returns true if x dominates y
function domine(x, y, opt::Optimisation=MAX)
    if opt == MIN
        return ((x[1] <= y[1] && x[2] < y[2])
            || (x[1] < y[1] && x[2] <= y[2])
            || (x[1] == y[1] && x[2] == y[2])) # No duplicates
    else
        return ((x[1] >= y[1] && x[2] > y[2])
            || (x[1] > y[1] && x[2] >= y[2])
            || (x[1] == y[1] && x[2] == y[2])) # No duplicates
    end
end

# ----- RATIOS --------------------------------------------------------------- #
# Computes the utilities (profit-to-weight ratios) for both objective functions
function utilities(prob::_MOMKP)

	p, n = size(prob.P)

	#=ratios = Matrix{Rational{Int}}(undef, p, n)
	for k in 1:p 
		for j in 1:n 
			ratios[k,j] = prob.P[k,j]//prob.W[1,j]
		end 
	end 
	return ratios=#

	r1 = [prob.P[1,j]//prob.W[1,j] for j in 1:n]
	r2 = [prob.P[2,j]//prob.W[1,j] for j in 1:n]
	return r1, r2
end

# ----- SOLUTIONS ------------------------------------------------------------ #
# Add an item to a solution
function addItem!(prob::_MOMKP, sol::Union{Solution,SolutionD}, item::Int)
	sol.X[item] = 1
	sol.z += prob.P[:,item]
end

# Add a break item to a solution
function addBreakItem!(prob::_MOMKP,
					   sol::Union{Solution,SolutionD},
					   ω_::Int,
					   item::Int)

	sol.X[item] = ω_//prob.W[1,item]
	sol.z += sol.X[item] * prob.P[:,item]
end

# Computes the dantzig solution for a given sequence
function dantzigSolution(prob::_MOMKP, sequence::Vector{Int})

	n   = size(prob.P)[2]
	ω_  = prob.ω[1]
	sol = Solution(n)
	i   = 1

	while i <= n && prob.W[1,sequence[i]] <= ω_
		item = sequence[i]
		# L'objet est inséré
		addItem!(prob, sol, item)
		ω_ -= prob.W[1,item]
		i += 1
	end

	return sol, i, ω_
end

# Builds a solution including the break item
function buildSolution(prob::_MOMKP, seq::Vector{Int})

	n          = size(prob.P)[2]
	sol, s, ω_ = dantzigSolution(prob, seq)

	if ω_ > 0
		# Une fraction de l'objet s est insérée
		addBreakItem!(prob, sol, ω_, seq[s])
	end

	return sol, s, ω_
end

# Re-build part of a solution after a sequence reversal
function reoptSolution(prob::_MOMKP,
					   seq::Vector{Int},
					   start::Int,
					   finish::Int, 
					   sol::Solution,
					   ω_::Int)

	# Les objets entre start et fin dans la séquence sont retirés
	for pos in start:finish

		item = seq[pos]
		
		if sol.X[item] < 1 && sol.X[item] > 0 
			# L'objet cassé est retiré
			sol.z      -= sol.X[item] * prob.P[:,item]
			sol.X[item] = 0
			
		elseif sol.X[item] == 1 
			# Un objet inséré dans le sac est retiré
			sol.z      -= prob.P[:,item]
			ω_         += prob.W[1,item]
			sol.X[item] = 0
		end 
	end

	# Les objets sont insérés en partant de start dans la nouvelle séquence
	pos = start
	while prob.W[1,seq[pos]] <= ω_

		# L'objet est inséré en entier
		addItem!(prob, sol, seq[pos])
		ω_ -= prob.W[1,seq[pos]]

		pos += 1
	end

	return sol, pos, ω_
end

# Returns true if the solution is integer 
function isInteger(sol::Solution, breakItem::Int)
	return (sol.X[breakItem] == 0 || sol.X[breakItem] == 1)
end

# ----- PROBLEMS ------------------------------------------------------------- #
# Groups together equivalent items 
function groupEquivalentItems(prob::_MOMKP)

    n      = size(prob.P)[2]
    r1, r2 = ratios(prob) 
    P1     = Int[] 
    P2     = Int[] 
    W      = Int[] 
    done   = Int[] # List of items that are equivalent to another item that has 
    # already been processed

    for i in 1:n       
        if !(i in done)

            p1 = prob.P[1,i]
            p2 = prob.P[2,i]
            w  = prob.W[1,i] 

            for j in i+1:n 

                if r1[i] == r1[j] && r2[i] == r2[j] 
                    # Items i and j are equivalent 
                    p1 += prob.P[1,j] 
                    p2 += prob.P[2,j] 
                    w  += prob.W[1,j] 
                    push!(done, j)
                end 
            end 

            push!(P1, p1)
            push!(P2, p2)
            push!(W, w) 
        end 
    end 

    P = Matrix(undef, 2, length(P1))
    for i in 1:length(P1) 
        P[1,i] = P1[i] ; P[2,i] = P2[i]  
    end 

    return _MOMKP(P, reshape(W, 1, length(W)), prob.ω)
end 

# ----- BOUND SETs ----------------------------------------------------------- #
# Returns true if sol is identical to the most recent solution in UB 
function identicalToPrevious(UB::DualBoundSet, sol::Solution)
    return length(UB.points) > 0 && UB.points[end] == sol.z 
end 

function identicalToPrevious(UB::DualBoundSet, y::Vector{Float64})
    return length(UB.points) > 0 && UB.points[end] == y
end 

# Updates the bound set by adding the point and corresponding constraint if it 
# is not already present and adding the solution to the list of integer 
# solutions if applicable 
function updateBoundSet!(upperBound::DualBoundSet, 
						 λ::Rational{Int}, 
						 sol::Solution, 
						 breakItem::Int)
    
    if !identicalToPrevious(upperBound, sol)
        push!(upperBound.points, sol.z) 
		push!(upperBound.constraints, Constraint(λ, sol.z))
        if isInteger(sol, breakItem) 
            push!(upperBound.integerSols, Solution(sol.X[1:end], sol.z[1:end])) 
        end 
    end 
end 

function updateBoundSet!(upperBound::DualBoundSet, 
						 λ::Union{Rational{Int},Float64}, 
						 y::Vector{Float64})

	if !identicalToPrevious(upperBound, y)
		push!(upperBound.points, y) 
		push!(upperBound.constraints, Constraint(λ, y))
	end 
end 

function updateBoundSet!(upperBound::DualBoundSet,  
						 sol::Solution, 
						 breakItem::Int)

	if !identicalToPrevious(upperBound, sol)
		push!(upperBound.points, sol.z) 
		if isInteger(sol, breakItem) 
			push!(upperBound.integerSols, Solution(sol.X[1:end], sol.z[1:end])) 
		end 
	end 
end 
