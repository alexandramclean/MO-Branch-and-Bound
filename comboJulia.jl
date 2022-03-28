
# Wrapper pour utiliser le code comboJulia.c/combo.c/combo.h dans Julia
# Première étape : compilation du code C (gcc -dynamiclib comboJulia.c -o comboJulia.dylib) dans un terminal standard

# Deuxième étape : fonction d'appel en Julia

struct item
    p::Int64 # profit
    w::Int64 # poids
    x::Int32 # {0,1} dans le sac ou pas
    i::Int32 # indice initial de l'objet (important! On appelle du code C, les indices commencent donc à 0)
end

comboJulia(prob::Vector{item},nbItem::Int32,capa::Int64,maxZ::Int64) =
ccall( (:comboJulia, "comboJulia"),
        Int64,
        (Ref{item}, Int32, Int64, Int64),
        prob, nbItem, capa, maxZ
        )

#= Sert plus à rien!
comboJuliaMoche(probp::Vector{Int64},probw::Vector{Int64},probx::Vector{Int32},probi::Vector{Int32},nbItem::Int32,capa::Int64,maxZ::Int64) =
ccall( (:comboJuliaMoche, "comboJulia"),
        Int64,
        (Ptr{Int64}, Ptr{Int64}, Ptr{Int32}, Ptr{Int32}, Int32, Int64, Int64),
        probp, probw, probx, probi, nbItem, capa, maxZ
        )
=#

# Essai avec le problème suivant
# Max z = 40x1 + 41x2 + 34x3 + 15x4 + 28x5
#   s.c.  16x1 + 12x2 + 11x3 + 5x4 + 16x5 <= 4
function essaiEvident1()
    prob::Vector{item} = Vector{item}(undef,5)
    prob[1] = item(40,16,0,0);
    prob[2] = item(41,12,0,1);
    prob[3] = item(34,11,0,2);
    prob[4] = item(15,5,0,3);
    prob[5] = item(28,16,0,4);
    capa::Int64 = 4
    maxZ::Int64 = 0
    for i in 1:5
        maxZ = maxZ + prob[i].p
    end
    nbItem::Int32 = 5
    val::Int32 = comboJulia(prob,nbItem,capa,maxZ)
    x::Vector{Int8} = Vector{Int8}(undef,5)
    index::Int64 = 0
    for i in 1:5
        index = prob[i].i + 1
        x[index] = prob[i].x
    end
    println("Valeur optimale = ",val)
    println("x = ",x)
end


# Essai avec le problème suivant
# Max z = 40x1 + 41x2 + 34x3 + 15x4 + 28x5
#   s.c.  16x1 + 12x2 + 11x3 + 5x4 + 16x5 <= 70
# Provoque un crash!!!
function essaiEvident2()
    prob::Vector{item} = Vector{item}(undef,5)
    prob[1] = item(40,16,0,0);
    prob[2] = item(41,12,0,1);
    prob[3] = item(34,11,0,2);
    prob[4] = item(15,5,0,3);
    prob[5] = item(28,16,0,4);
    capa::Int64 = 70
    maxZ::Int64 = 0
    for i in 1:5
        maxZ = maxZ + prob[i].p
    end
    nbItem::Int32 = 5
    val::Int32 = comboJulia(prob,nbItem,capa,maxZ)
    x::Vector{Int8} = Vector{Int8}(undef,5)
    index::Int64 = 0
    for i in 1:5
        index = prob[i].i + 1
        x[index] = prob[i].x
    end
    println("Valeur optimale = ",val)
    println("x = ",x)
end


# Essai avec le problème suivant
# Max z = 40x1 + 41x2 + 34x3 + 15x4 + 28x5
#   s.c.  16x1 + 12x2 + 11x3 + 5x4 + 16x5 <= 30
#         xj in {0,1}, j = 1..5

function essai()
    prob::Vector{item} = Vector{item}(undef,5)
    prob[1] = item(40,16,0,0);
    prob[2] = item(41,12,0,1);
    prob[3] = item(34,11,0,2);
    prob[4] = item(15,5,0,3);
    prob[5] = item(28,16,0,4);

    capa::Int64 = 30
    maxZ::Int64 = 0
    for i in 1:5
        maxZ = maxZ + prob[i].p
    end

    nbItem::Int32 = 5
    val::Int32 = comboJulia(prob,nbItem,capa,maxZ)
    x::Vector{Int8} = Vector{Int8}(undef,5)
    index::Int64 = 0
    for i in 1:5
        index = prob[i].i + 1
        x[index] = prob[i].x
    end
    println("Valeur optimale = ",val)
    println("x = ",x)
end
