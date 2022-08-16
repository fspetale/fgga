fgga <- function(graphOnto, tableOntoTerms, dxCharacterized,
                dxTestCharacterized, kFold, kernelSVM = "radial", tmax = 200,
                epsilon = 0.001) {
    ontoNodes <- nodes(graphOnto)
    rootOnto <- leaves(graphOnto, "in")
    modelSVMs <- lapply(ontoNodes,
        FUN = svmTrain, tableOntoTerms = tableOntoTerms,
        dxCharacterized = dxCharacterized, graphOnto = graphOnto,
        kernelSVM = kernelSVM)

    varianceOntoTerms <- varianceOnto(tableOntoTerms, dxCharacterized,
                                kFold, graphOnto, rootOnto, kernelSVM)

    matrixOntoTest <- svmOnto(svmMoldel = modelSVMs,
                                dxCharacterized = dxTestCharacterized,
                                rootNode = rootOnto,
                                varianceSVM = varianceOntoTerms)

    modelFG <- fgga2bipartite(graphOnto)

    matrixFGGATest <- t(apply(matrixOntoTest, MARGIN = 1, FUN = msgFGGA,
        matrixFGGA = modelFG, graphOnto = graphOnto, tmax = tmax,
        epsilon = epsilon))

    return(matrixFGGATest)
}
