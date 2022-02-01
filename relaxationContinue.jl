################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Calcul paramétrique de la relaxation continue                        #
################################################################################

include("parserMomkpPG.jl")
include("parserMomkpZL.jl")

mutable struct Solution
    X::Vector{Union{Int64,Float64}} # Liste des éléments insérés dans le sac
    z::Vector{Float64}              # Valeurs pour les fonctions objectif
end
Solution(n) = Solution(zeros(n), [0,0])

function copy(sol::Solution)
    return Solution(sol.X[1:end],
                    sol.z[1:end])
end

# Précondition : cas uni-dimensionnel
function ratios(prob::_MOMKP)

	n  = size(prob.P)[2]
	r1 = Vector{Float64}(undef, n)
	r2 = Vector{Float64}(undef, n)
	for i in 1:n 
		r1[i] = prob.P[1,i]/prob.W[1,i] 
		r2[i] = prob.P[2,i]/prob.W[1,i]
	end
	return r1, r2
end

# Calcul des poids critiques
function poidsCritiques(prob::_MOMKP, r1, r2)

	n       = size(prob.P)[2]
	poids   = Float64[]
	couples = Tuple{Int64,Int64}[] 
	# Calcul des poids critiques pour chaque couple d'objets (i,j)
	for i in 1:n
		for j in i:n
			λ = (r2[j] - r2[i])/(r1[i]-r2[i]-r1[j]+r2[j])
			if λ > 0 && λ < 1
				push!(poids, λ)
				push!(couples, (i,j))
			end
		end
	end
	# Tri des poids critiques dans l'ordre décroissant
	perm = sortperm(poids, rev=true)
	return poids[perm], couples[perm]
end

function consSolution(prob::_MOMKP, sequence)

	n              = size(prob.P)[2]
	capaResiduelle = prob.ω[1]
	sol            = Solution(n)
	i              = 1
	
	while capaResiduelle > 0 
		objet = sequence[i]
		
		if prob.W[1,objet] <= capaResiduelle # L'objet est inséré en entier
			sol.X[objet] = 1
			sol.z += prob.P[:,objet]
			capaResiduelle -= prob.W[1,objet]
			
		else # Une fraction de l'objet est insérée
			sol.X[objet] = capaResiduelle/prob.W[1,objet] 
			sol.z += sol.X[objet]*prob.P[:,objet]
		end
		i += 1
	end
	
	# Position de l'objet cassé 
	if sol.X[sequence[i-1]] < 1
		s = i-1 
	else 
		s = i
	end
	return sol, s, capaResiduelle
end

function transposition(prob::_MOMKP, seq, sol_init, capaResiduelle, i, j)	
	
	sol = copy(sol_init)
	
	# L'objet à la position i est enlevé du sac
	sol.X[seq[i]] = 0
	sol.z -= prob.P[:,seq[i]]
	
	# La fraction de l'objet à la position j est enlevée
	if sol.X[seq[j]] < 1 && sol.X[seq[j]] > 0
		sol.z -= sol.X[seq[j]]*prob.P[:,seq[j]]
	end
			
	if prob.W[1,seq[j]] <= capaResiduelle
	# L'objet à la position j peut être entièrement inséré
		sol.X[seq[j]] = 1
		sol.z += prob.P[:,seq[j]]
		capaResiduelle -= prob.W[1,seq[j]] 
				
		if capaResiduelle > 0 # Il reste de la place dans le sac
			# L'objet à la place i devient l'objet cassé 
			sol.X[seq[i]] = capaResiduelle/prob.W[1,seq[i]] 
			sol.z += sol.X[seq[i]]*prob.P[:,seq[i]] 
		end 
		
	else # L'objet j est cassé 
		sol.X[seq[j]] = capaResiduelle/prob.W[1,seq[j]] 
		sol.z += sol.X[seq[j]]*prob.P[:,seq[j]] 
	end 
	return sol	
end 		

function relaxationContinue(prob::_MOMKP)

	solutionsObtenues = Solution[]

	# Calcul des ratios
	r1, r2 = ratios(prob)
	
	# Calcul des poids critiques
	poids, couples = poidsCritiques(prob, r1, r2)
	
	# Tri des ratios dans l'ordre lexicographique décroissant selon (r1,r2)
	seq = sortperm(r1, rev=true) # Séquence d'objets
	pos = sortperm(seq)     # Position des objets dans la séquence
	
	println("Séquence : ", seq)
	println("Positions : ", pos)
	
	# Construction de la première solution
	sol, s, capaResiduelle = consSolution(prob, seq)
	push!(solutionsObtenues, sol)
	
	# Boucle principale
	for iter in 1:length(poids)
	
		println("\nIter ", iter)
		
		(i,j) = couples[iter] 
		k = min(pos[i], pos[j])

		# La place de l'objet cassé est échangée avec un objet dans le sac
		# OU La place d'un objet dans le sac est échangée avec un objet qui n'est pas 
		# dans le sac
		if k == s-1 
			capaResiduelle += prob.W[1,seq[k]] # L'objet à la position k est retiré
			sol = transposition(prob, seq, sol, capaResiduelle, s-1, s)
		
			if sol.X[seq[k+1]] < 1 && sol.X[seq[k+1]] > 0
				println("Ici")
				s = s-1
			end 
		# La place de l'objet cassé est échangée avec un objet qui n'est pas 
		# dans le sac
		elseif k == s 
			sol = transposition(prob, seq, sol, capaResiduelle, s, s+1)
			if sol.X[seq[k+1]] == 1
				println("Là")
				s = s+1
			end
		end 
			
		# Mise à jour de la séquence
		tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
		seq[pos[i]] = i ; seq[pos[j]] = j
		
		println("Séquence : ", seq)
		println("Position de l'objet cassé : ", s)
		println("Solution : ", sol)
		
		push!(solutionsObtenues, sol) 
	end
	
	return solutionsObtenues		

end

# Exemple 
didactic = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])
solutionsObtenues = relaxationContinue(didactic)

for sol in solutionsObtenues 
	println(sol.z)
end
