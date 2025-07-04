################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Parsers                                                              #
################################################################################

include("parserMomkpPG.jl")
include("parserMomkpZL.jl")
include("newDPKP.jl")

# Function for parsing a file containing a reference set
function readReferenceSet(fname::String)

    f = open(fname)
    lines = readlines(f)

    # Reference set 
    ref = Vector{Vector{Float64}}(undef, length(lines)-2)

    for i in 3:length(lines)
        splitLine = split(lines[i])
        ref[i-2] = [parse(Float64, splitLine[1]), parse(Float64, splitLine[2])]
    end 

    return ref
end 

function readInstanceKP(fname::String)

    prob = parseKP(fname)
    return _MOMKP(transpose(hcat(prob.p1,prob.p2)), 
                  reshape(prob.w, 1, :),
                  [prob.Omega])
end 