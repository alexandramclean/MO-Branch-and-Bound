################################################################################
#         Stage M2 - Ensembles bornants pour un B&B multi-objectif             #
#         MCLEAN Alexandra, encadré par Gandibleux X. et Przybylski A.         #
#         Branch-and-bound algorithm                                           #
################################################################################

include("branchAndBoundFunctions.jl")

# ----- BRANCH-AND-BOUND ----------------------------------------------------- #
# Recursive branching function 
function branch!()
end 

function branchAndBound(prob::_MOMKP) 

    rootNode = Node(DualBoundSet(), Nothing, NOTPRUNED)

end 