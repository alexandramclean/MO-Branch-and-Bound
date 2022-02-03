# oh non
include("functions.jl")
include("parserMomkpPG.jl")
include("parserMomkpZL.jl")

function transpositions(prob::_MOMKP, verbose=false)

	# Calcul des ratios
	r1, r2 = ratios(prob)
	
	# Calcul des poids critiques
	weights, pairs = criticalWeights(prob, r1, r2)
	
	# Tri des ratios dans l'ordre lexicographique décroissant selon (r1,r2)
	seq = sortperm(r1, rev=true) # Séquence d'objets
	pos = sortperm(seq)          # Position des objets dans la séquence
	
	if verbose 
		println(seq) 
	end
	
	# Boucle principale
	for iter in 1:length(weights)
		
		#println("\nIter ", iter)
		
		(i,j) = pairs[iter] 
		k = min(pos[i], pos[j])
		
		if iter < length(weights) && weights[iter] == weights[iter+1] 
			println("égalité")
		end 
		
		if !(pos[i] == pos[j]+1 || pos[j] == pos[i]+1) 
			println("oh non") 
		end 
		
		# Mise à jour de la séquence
		tmp = pos[i] ; pos[i] = pos[j] ; pos[j] = tmp
		seq[pos[i]] = i ; seq[pos[j]] = j
		
		if verbose 
			println(seq)
		end
	end
end

#prob = _MOMKP([11 2 8 10 9 1 ; 2 7 8 4 1 3], [4 4 6 4 3 2], [11])

# Exemple d'instances où l'on supprime beaucoup de paires d'objets 
#prob = _MOMKP([6 2 4 1 7 ; 2 4 1 3 5], [1 3 2 4 4], [12])


fname = "../instancesPG/set2/A2.DAT"
if fname[length(fname)-3:length(fname)] == ".DAT"
    prob = readInstanceMOMKPformatPG(false, fname)
else
    prob = readInstanceMOMKPformatZL(false, fname)
end

r1, r2 = ratios(prob)

weights, pairs = criticalWeights(prob, r1, r2)
println("pairs : ", pairs, "\n")

transpositions(prob) 

