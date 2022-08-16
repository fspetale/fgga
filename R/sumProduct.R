.productMsg <- function(valuesTPG, matrixFunction, matrixNode, index1, index2) {
    valueProduct <- 1
    for (z in seq_len(length(valuesTPG) - 1)) {
        valuesMsg <- matrixNode[matrixFunction[index1, index2][[1]]$neighbor[z],
                    index1][[1]][[valuesTPG[z]]]
        valueProduct <- valuesMsg * valueProduct
    }
    return(valueProduct)
}

msgFGGA <- function(matrixFGGA, obsValueOntoTerms, graphOnto, tmax = 200,
                    epsilon = 0.001) {
    if (!is(matrixFGGA, 'fgga')){
        message("matrixFGGA is not class fgga")
        return()
    } else {
        matrixFGGA <- as(matrixFGGA, "matrix")
    }
    nFactors <- dim(matrixFGGA)[2]
    nNodes <- dim(matrixFGGA)[1] / 2
    matrixOnto <- as(graphOnto, "matrix")
    matrixNodeToFunction <- array(data = list(), dim = c(nNodes, nFactors))
    matrixFunctionToNode <- array(data = list(), dim = c(nFactors, nNodes))
    rownames(matrixNodeToFunction) <- rownames(matrixFGGA)[seq_len(nNodes)]
    colnames(matrixFunctionToNode) <- rownames(matrixFGGA)[seq_len(nNodes)]
    colnames(matrixNodeToFunction) <- colnames(matrixFGGA)
    rownames(matrixFunctionToNode) <- colnames(matrixFGGA)

    # Initialization
    for (i in seq_len(dim(matrixFunctionToNode)[2])) {
        matrixFunctionToNode[nNodes - 1 + i, i][[1]] <- list(
            zeroValues = 1 - obsValueOntoTerms[i],
            oneValues = obsValueOntoTerms[i], valueOnto = 1)
        if (i < nNodes) {
            nodeLinks <- which(matrixFGGA[, i] == 1)
            submatrixOnto <- matrixOnto[nodeLinks, nodeLinks]
            leafOnto <- names(which(apply(submatrixOnto, MARGIN = 1, sum) == 0))
            parentsOnto <- setdiff(colnames(submatrixOnto), leafOnto)
            functionLink <- t(apply(t(names(nodeLinks)), MARGIN=2,
                FUN = function(x, y, w) {(if (x == w) setdiff(y, x) else c(
                w, setdiff(y, c(x, w))))}, w = leafOnto, y = names(nodeLinks)))
            for (k in seq_len(length(nodeLinks))) {
                    if (length(nodeLinks) > 2) {
                        matrixFunctionToNode[i, nodeLinks[k]][[1]] <-
                            list(zeroValues = 1, oneValues = 1,
                            neighbor = functionLink[k, ], leafOnto = leafOnto,
                            parentsOnto = parentsOnto,
                            TPG = length(nodeLinks) - 1)
                    } else {
                        matrixFunctionToNode[i, nodeLinks[k]][[1]] <- list(
                            zeroValues = 1,
                            oneValues = 1, neighbor = functionLink[1, k],
                            leafOnto = leafOnto, parentsOnto = parentsOnto,
                            TPG = length(nodeLinks) - 1)
                    }
            }
        }
        nodeLinks <- which(matrixFGGA[i, ] == 1)
        functionLink <- t(apply(t(names(nodeLinks)), MARGIN=2,
            FUN = function(x, y) {(setdiff(y, x))}, y = names(nodeLinks)))
        for (k in seq_len(length(nodeLinks))) {
                if (length(nodeLinks) > 2) {
                    matrixNodeToFunction[i, nodeLinks[k]][[1]] <- list(
                    zeroValues = 1, oneValues = 1, neighbor = functionLink[k, ])
                } else {
                    matrixNodeToFunction[i, nodeLinks[k]][[1]] <- list(
                    zeroValues = 1, oneValues = 1, neighbor = functionLink[1, k]
                    )
                }
        }
    }
    rm(submatrixOnto, leafOnto, parentsOnto)
    for (t in seq_len(tmax)) {
        # Msg node variable to function
        for (i in seq_len(nNodes)) {
            activeFunction <- which(matrixFGGA[i, ] == 1)
            for (j in seq_len(length(activeFunction))) {
                valuesMsg <- apply(t(activeFunction[-j]), MARGIN = 2,
                    FUN = function(x, y, z) (y[x, z][[1]]$oneValues),
                    y = matrixFunctionToNode, z = i)
                matrixNodeToFunction[i, activeFunction[j]][[1]]$oneValues <-
                    prod(valuesMsg)
                valuesMsg <- apply(t(activeFunction[-j]), MARGIN = 2,
                    FUN = function(x, y, z) (y[x, z][[1]]$zeroValues),
                    y = matrixFunctionToNode, z = i)
                matrixNodeToFunction[i, activeFunction[j]][[1]]$zeroValues <-
                    prod(valuesMsg)
                gammaMsg <- 1 / (matrixNodeToFunction[i,
                activeFunction[j]][[1]]$oneValues + matrixNodeToFunction[i,
                        activeFunction[j]][[1]]$zeroValues)
                matrixNodeToFunction[i, activeFunction[j]][[1]]$oneValues <-
                    matrixNodeToFunction[i,
                    activeFunction[j]][[1]]$oneValues * gammaMsg
                matrixNodeToFunction[i, activeFunction[j]][[1]]$zeroValues <-
                    matrixNodeToFunction[i,
                    activeFunction[j]][[1]]$zeroValues * gammaMsg
            }
        }
        # Msg function to node variable
        for (i in seq_len(nNodes - 1)) {
            activeFunction <- which(matrixFGGA[, i] == 1)
            for (j in seq_len(length(activeFunction))) {
                functionsTPG <- tableTPG(matrixFunctionToNode[i,
                                        activeFunction[j]][[1]]$TPG+1)
                if (matrixFunctionToNode[i, activeFunction[j]][[1]]$leafOnto ==
                    colnames(matrixFunctionToNode)[activeFunction[j]]) {
                    for (k in seq(1, 2)) {
                        indexMsg <- which(functionsTPG[, 1] == k)
                        subMatrixTPG <- functionsTPG[indexMsg, -1]
                        if (is.matrix(subMatrixTPG)) {
                            valuesMsg <- apply(subMatrixTPG,
                                MARGIN = 1, FUN = .productMsg,
                                index1 = i, index2 = activeFunction[j],
                                matrixFunction = matrixFunctionToNode,
                                matrixNode = matrixNodeToFunction)
                        } else {
                                valuesMsg <- apply(as.matrix(t(subMatrixTPG)),
                                    MARGIN = 1, FUN = .productMsg, index1 = i,
                                    index2 = activeFunction[j],
                                    matrixFunction = matrixFunctionToNode,
                                    matrixNode = matrixNodeToFunction)
                            }
                        if (k == 1) {
                            matrixFunctionToNode[i,
                                activeFunction[j]][[1]]$zeroValues <-
                                sum(valuesMsg)
                        } else {
                            matrixFunctionToNode[i,
                                activeFunction[j]][[1]]$oneValues <-
                                sum(valuesMsg)
                        }
                    }
                } else {
                    for (k in seq(1, 2)) {
                        indexMsg <- which(functionsTPG[, 2] == k)
                        subMatrixTPG <- functionsTPG[indexMsg, -2]
                        if (is.matrix(subMatrixTPG)) {
                            valuesMsg <- apply(subMatrixTPG, MARGIN = 1,
                                    FUN = .productMsg, index1 = i,
                                    index2 = activeFunction[j],
                                    matrixFunction = matrixFunctionToNode,
                                    matrixNode = matrixNodeToFunction)
                        } else {
                                valuesMsg <- apply(as.matrix(t(subMatrixTPG)),
                                MARGIN = 1, FUN = .productMsg, index1 = i,
                                index2 = activeFunction[j],
                                matrixFunction = matrixFunctionToNode,
                                matrixNode = matrixNodeToFunction)
                        }
                        if (k == 1) {
                            matrixFunctionToNode[i,
                            activeFunction[j]][[1]]$zeroValues <- sum(valuesMsg)
                        } else {
                            matrixFunctionToNode[i,
                            activeFunction[j]][[1]]$oneValues <- sum(valuesMsg)}
                    }
                }
                gammaMsg <- 1 / (matrixFunctionToNode[i,
                    activeFunction[j]][[1]]$oneValues +
                    matrixFunctionToNode[i, activeFunction[j]][[1]]$zeroValues)
                matrixFunctionToNode[i, activeFunction[j]][[1]]$oneValues <-
                    matrixFunctionToNode[i, activeFunction[j]][[1]]$oneValues *
                        gammaMsg
                matrixFunctionToNode[i, activeFunction[j]][[1]]$zeroValues <-
                    matrixFunctionToNode[i, activeFunction[j]][[1]]$zeroValues *
                        gammaMsg
            }
        }
        ValuesOnto <- lapply(diag(matrixNodeToFunction[,
            seq(nNodes, nFactors)]), FUN = function(x) (x[c(1, 2)]))
        oneValuesOnto <- vapply(seq(1, nNodes), FUN = function(x, y, z) {
            (y[[x]]$oneValues * z[x])}, y = ValuesOnto, z = obsValueOntoTerms,
            FUN.VALUE = numeric(1))
        zeroValuesOnto <- vapply(seq(1, nNodes), FUN = function(x, y, z) {
            (y[[x]]$zeroValues * (1 - z[x]))}, y = ValuesOnto,
            z = obsValueOntoTerms, FUN.VALUE = numeric(1))
        oneValuesOnto <- vapply(seq(1, nNodes), FUN = function(x, y, z) {
            (y[x] / (y[x] + z[x]))}, y = oneValuesOnto, z = zeroValuesOnto,
            FUN.VALUE = numeric(1))
        ValuesOnto <- vapply(diag(matrixFunctionToNode[
            seq(nNodes, nFactors), ]), FUN = function(x) (x$valueOnto),
            FUN.VALUE = numeric(1))
        errorGO <- vapply(seq(1, nNodes), FUN = function(x, y, z) {
            (abs(y[x] - z[x]))}, y = oneValuesOnto, z = ValuesOnto,
            FUN.VALUE = numeric(1))
        if (max(errorGO) < epsilon) {
            names(oneValuesOnto) <- colnames(matrixFunctionToNode)
            return(oneValuesOnto)
        } else {
            for (i in seq_len(nNodes)) {
                matrixFunctionToNode[nNodes - 1 + i, i][[1]]$valueOnto <-
                    oneValuesOnto[i]}
        }
    }
    ValuesOnto <- array(0, nNodes)
    names(ValuesOnto) <- colnames(matrixFunctionToNode)
    message("Graph UNCONVERGED")
    return(ValuesOnto)
}
