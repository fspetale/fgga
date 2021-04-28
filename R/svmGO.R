svmTrain <- function(nodeGraph, tableGOs, dxCharacterized, graphGO,
                        kernelSVM = "radial") {
    print(nodeGraph)
    tableGOs <- tableGOs[order(rownames(tableGOs)), ]
    dxCharacterized <- dxCharacterized[order(rownames(dxCharacterized)), ]
    modelSvm_reduce <- c()
    zeroValue <- array(0, 0.1 * dim(dxCharacterized)[2])
    idAncestor <- ancestors(nodeGraph, graphGO)
    rootNodes <- c("GO:0008150", "GO:0003674", "GO:0005575", "GO:FES")
    rootNodes <- intersect(idAncestor, rootNodes)
    idAncestor <- setdiff(idAncestor, rootNodes)
    indices <- which(tableGOs[, nodeGraph] == 1)
    positiveGeneNames <- rownames(tableGOs)[indices]
    notNegativeGeneNames <- positiveGeneNames
    if (length(idAncestor) > 0) {
        for (i in seq_len(length(idAncestor))) {
            buffer <- rownames(tableGOs)[which(tableGOs[, idAncestor[i]] == 1)]
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
                id.negative <- sample(seq_len(length(negativeGeneNames)), nn)
                exprsValues <- dxCharacterized[c(
                    positiveGeneNames,
                    negativeGeneNames[id.negative]
                ), ]
            }
            else {
                np <- 2 * nn
                id.positive <- sample(seq_len(length(positiveGeneNames)), np)
                exprsValues <- dxCharacterized[c(
                    positiveGeneNames[id.positive],
                    negativeGeneNames
                ), ]
            }
        }
        else {
            exprsValues <- dxCharacterized[c(
                positiveGeneNames,
                negativeGeneNames
            ), ]
        }
        exprsValues[is.na(exprsValues)] <- 0
        gos <- as.factor(c(rep(1, np), rep(2, nn)))
        zeroValue <- apply(exprsValues, MARGIN = 2, mean)
        zeroValue <- which(zeroValue == 0)
        if (length(zeroValue) != 0) exprsValues <- exprsValues[, -zeroValue]
        modelSvm <- svm(exprsValues, gos,
            method = "C-classification", cost = 500,
            kernel = kernelSVM, na.action = na.omit
        )
        modelSvm_reduce <- modelSvm[-c(28, 29)]
        class(modelSvm_reduce) <- class(modelSvm)
        modelSvm_reduce$zeroValue <- zeroValue
    }
    modelSvm_reduce$goTerm <- nodeGraph
    return(modelSvm_reduce)
}

svmGO <- function(svmMoldel, dxCharacterized, rootNode, varianceSVM) {
    dxTestSVM <- lapply(svmMoldel,
        FUN = .svmTest, dxTest = dxCharacterized,
        rootGO = rootNode, variance = varianceSVM
    )
    matrixTest <- vapply(dxTestSVM, cbind,
                        FUN.VALUE = double(dim(dxTestSVM[[1]])[1]))
    colnames(matrixTest) <- vapply(dxTestSVM, colnames,
                                   FUN.VALUE = character(1))
    rownames(matrixTest) <- vapply(dxTestSVM[1], rownames,
                                FUN.VALUE = character(dim(dxTestSVM[[1]])[1]))
    return(matrixTest)
}

.svmTest <- function(svmMoldel, dxTest, rootGO, variance) {
    dxTest <- as.data.frame(dxTest)
    matrixScore <- matrix(0, dim(dxTest)[1], 1)
    rownames(matrixScore) <- rownames(dxTest)
    nodeGO <- svmMoldel$goTerm
    colnames(matrixScore) <- nodeGO
    if (nodeGO != rootGO) {
        zeroValues <- svmMoldel$zeroValue
        if (length(zeroValues) != 0) dxTest <- dxTest[, -zeroValues]
        predGO <- predict(svmMoldel, dxTest, decision.values = TRUE)
        matrixScore[, 1] <- 1 / (1 + exp(-2 * attr(predGO, "decision.values") /
            variance[nodeGO]))
    } else {
        matrixScore[, 1] <- 0.9999
    }
    matrixScore[which(matrixScore < 1.0e-06)] <- 1.0e-06
    matrixScore[which(matrixScore > 0.9999999)] <- 0.9999
    return(matrixScore)
}
