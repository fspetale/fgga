varianceGO <- function(tableGOs, dxCharacterized, kFold, graphGO, rootNode,
                        kernelSVM = "radial") {
    tableGOs <- as.matrix(tableGOs[order(rownames(tableGOs)), ])
    dxCharacterized <- as.data.frame(dxCharacterized[order(rownames(
        dxCharacterized
    )), ])
    nodesGraph <- colnames(tableGOs)
    variance <- array(0, length(nodesGraph))
    names(variance) <- nodesGraph
    for (j in seq_len(length(nodesGraph))) {
        print(nodesGraph[j])
        if (nodesGraph[j] != rootNode) {
            checkPosiSample <- length(which(tableGOs[, j] == 1))
            if (checkPosiSample > 1) {
                if (checkPosiSample > kFold) {
                    nFold <- kFold
                } else {
                    nFold <- checkPosiSample
                }
                dxFold <- createFolds(as.factor(tableGOs[, j]), nFold)
                varianceList <- c()
                for (i in seq_len(kFold)) {
                    varianceList[i] <- .varianceSVM(
                        indexFold = i, dxFolds = dxFold, tableGOs = tableGOs,
                        dxCharacterized = dxCharacterized,
                        nodeGO = nodesGraph[j],
                        graphGO = graphGO, kernelSVM = kernelSVM)
                }
                variance[j] <- sum(varianceList)
            } else {
                variance[j] <- 1
            }
            variance[j] <- variance[j] / (nFold)
        }
    }
    variance[which(variance < 1e-06)] <- 0.00001
    return(variance)
}

.varianceSVM <- function(indexFold, dxFolds, tableGOs, dxCharacterized, nodeGO,
                            graphGO, kernelSVM = "radial") {
    targetGO <- tableGOs[dxFolds[[indexFold]], nodeGO]
    dxCharacterizedValid <- dxCharacterized[dxFolds[[indexFold]], ]
    tableGOs <- tableGOs[-dxFolds[[indexFold]], ]
    dxCharacterized <- dxCharacterized[-dxFolds[[indexFold]], ]
    zeroValues <- array(0, 0.1 * dim(dxCharacterized)[2])
    ancestorID <- ancestors(nodeGO, graphGO)
    rootNodes <- c("GO:0008150", "GO:0003674", "GO:0005575", "GO:FES")
    rootNodes <- intersect(ancestorID, rootNodes)
    ancestorID <- setdiff(ancestorID, rootNodes)
    indices <- which(tableGOs[, nodeGO] == 1)
    positiveGeneNames <- rownames(tableGOs)[indices]
    notNegativeGeneNames <- positiveGeneNames
    if (length(ancestorID) > 0) {
        for (i in seq_len(length(ancestorID))) {
            buffer <- rownames(tableGOs)[which(tableGOs[, ancestorID[i]] == 1)]
            notNegativeGeneNames <- union(notNegativeGeneNames, buffer)
        }
    }

    negativeGeneNames <- setdiff(rownames(tableGOs), notNegativeGeneNames)
    np <- length(positiveGeneNames)
    nn <- length(negativeGeneNames)
    if ((np != 0) && (nn != 0)) {
        if (nn > 2 * np || np > 2 * nn) {
            if (nn > 2 * np) {
                nn <- 2 * np
                idNegative <- sample(seq_len(length(negativeGeneNames)), nn)
                exprsValues <- dxCharacterized[
                    c(positiveGeneNames, negativeGeneNames[idNegative]), ]
            } else {
                np <- 2 * nn
                idPositive <- sample(seq_len(length(positiveGeneNames)), np)
                exprsValues <- dxCharacterized[
                    c(positiveGeneNames[idPositive], negativeGeneNames), ]
            }
        } else {
            exprsValues <- dxCharacterized[c(
                positiveGeneNames, negativeGeneNames), ]
        }
        exprsValues[is.na(exprsValues)] <- 0
        gos <- as.factor(c(rep(1, np), rep(2, nn)))
        zeroValues <- apply(exprsValues, MARGIN = 2, mean)
        zeroValues <- which(zeroValues == 0)
        if (length(zeroValues) != 0) exprsValues <- exprsValues[, -zeroValues]

        model_svm <- svm(exprsValues, gos,
            method = "C-classification", cost = 500,
            kernel = kernelSVM, na.action = na.omit)

        if (length(zeroValues) != 0) {
                dxCharacterizedValid <- dxCharacterizedValid[, -zeroValues]
            }
        predGO <- predict(model_svm, dxCharacterizedValid,
                            decision.values = TRUE)
        targetGO[which(targetGO == 0)] <- (-1)
        w <- which(targetGO == 1)
        error <- sum((attr(predGO, "decision.values")[w] -
            mean(attr(predGO, "decision.values")[w]))^2)
        variancePos <- error / length(w)
        error <- sum((attr(predGO, "decision.values")[-w] -
            mean(attr(predGO, "decision.values")[-w]))^2)
        varianceNeg <- error / (length(predGO) - length(w))
    } else {
        variancePos <- 0
    }
    return(variancePos)
}
