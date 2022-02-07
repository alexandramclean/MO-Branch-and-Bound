include("functions.jl")

fname = "../instancesPG/set2/A2.DAT"

if fname[length(fname)-3:length(fname)] == ".DAT"
    prob = readInstanceMOMKPformatPG(false, fname)
else
    prob = readInstanceMOMKPformatZL(false, fname)
end

r1, r2 = ratios(prob)
n = length(r1)

for i in 1:n 
	for j in i:n 
		numerateur, denominateur = lambda(r1, r2, i, j) 
		if denominateur != 0 
			λ = numerateur//denominateur 
			println("(", i, ",", j, ") : ", "λ = ", λ, " ~ ", numerateur/denominateur*1.0)
		end
	end
end
