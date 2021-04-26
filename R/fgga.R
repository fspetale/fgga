fgga <- function(graphGO, tableGOs, dxCharacterized, dxTestCharacterized,
                 kFold, kernelSVM = "radial", tmax = 200, epsilon = 0.001){

    goNodes <- nodes(graphGO)
    rootGO <- leaves(graphGO, "in")
    modelSVMs <- lapply(goNodes, FUN = svmTrain, tableGOs = tableGOs,
                        dxCharacterized = dxCharacterized,
                        graphGO = graphGO, kernelSVM = kernelSVM)

    varianceGOs <- varianceGO(tableGOs, dxCharacterized, kFold, graphGO,
                              rootGO, kernelSVM)

    matrixGOTest <- svmGO(svmMoldel = modelSVMs,
                          dxCharacterized = dxTestCharacterized,
                          rootNode = rootGO, varianceSVM = varianceGOs)

    modelFG <- fgga2bipartite(graphGO)

    matrixFGGATest <- t(apply(matrixGOTest, MARGIN = 1, FUN = msgFGGA,
                              matrixFGGA = modelFG, graphGO = graphGO,
                              tmax = tmax, epsilon = epsilon))

    return(matrixFGGATest)
    }
