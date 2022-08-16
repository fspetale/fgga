svmTrain <- function(nodeGraph, tableOntoTerms, dxCharacterized, graphOnto,
                        kernelSVM = "radial") {
    message("Term: ", nodeGraph)
    tableOntoTerms <- tableOntoTerms[order(rownames(tableOntoTerms)), ]
    dxCharacterized <- dxCharacterized[order(rownames(dxCharacterized)), ]
    modelSvm_reduce <- c()
    zeroValue <- array(0, 0.1 * dim(dxCharacterized)[2])
    idAncestor <- ancestors(nodeGraph, graphOnto)
    rootNodes <- c("GO:0008150", "GO:0003674", "GO:0005575", "PO:0025131",
                    "PO:0009012", "ZFA:0100000", "HP:0000001", "Ontology:FES")
    rootNodes <- intersect(idAncestor, rootNodes)
    idAncestor <- setdiff(idAncestor, rootNodes)
    indices <- which(tableOntoTerms[, nodeGraph] == 1)
    positiveGeneNames <- rownames(tableOntoTerms)[indices]
    notNegativeGeneNames <- positiveGeneNames
    if (length(idAncestor) > 0) {
        for (i in seq_len(length(idAncestor))) {
            buffer <- rownames(tableOntoTerms)[which(tableOntoTerms[,
                        idAncestor[i]] == 1)]
            notNegativeGeneNames <- union(notNegativeGeneNames, buffer)
        }
    }

    negativeGeneNames <- setdiff(rownames(tableOntoTerms), notNegativeGeneNames)
    np <- length(positiveGeneNames)
    nn <- length(negativeGeneNames)
    if ((np != 0) && (nn != 0)) {
        if (nn > 2 * np || np > 2 * nn) {
            if (nn > 2 * np) {
                nn <- 2 * np
                id.negative <- sample(seq_len(length(negativeGeneNames)), nn)
                exprsValues <- dxCharacterized[c(
                    positiveGeneNames, negativeGeneNames[id.negative]), ]
            }
            else {
                np <- 2 * nn
                id.positive <- sample(seq_len(length(positiveGeneNames)), np)
                exprsValues <- dxCharacterized[c(
                    positiveGeneNames[id.positive], negativeGeneNames), ]
            }
        }
        else {
            exprsValues <- dxCharacterized[c(
                positiveGeneNames, negativeGeneNames), ]
        }
        exprsValues[is.na(exprsValues)] <- 0
        terms <- as.factor(c(rep(1, np), rep(2, nn)))
        zeroValue <- apply(exprsValues, MARGIN = 2, mean)
        zeroValue <- which(zeroValue == 0)
        if (length(zeroValue) != 0) exprsValues <- exprsValues[, -zeroValue]
        modelSvm <- svm(exprsValues, terms, method = "C-classification",
                        cost = 500, kernel = kernelSVM, na.action = na.omit)
        modelSvm_reduce <- modelSvm[-c(28, 29)]
        class(modelSvm_reduce) <- class(modelSvm)
        modelSvm_reduce$zeroValue <- zeroValue
    }
    modelSvm_reduce$ontoTerm <- nodeGraph
    return(modelSvm_reduce)
}

svmOnto <- function(svmMoldel, dxCharacterized, rootNode, varianceSVM) {
    dxTestSVM <- lapply(svmMoldel, FUN = .svmTest, dxTest = dxCharacterized,
        rootOnto = rootNode, variance = varianceSVM)
    matrixTest <- vapply(dxTestSVM, cbind,
                        FUN.VALUE = double(dim(dxTestSVM[[1]])[1]))
    colnames(matrixTest) <- vapply(dxTestSVM, colnames,
                                    FUN.VALUE = character(1))
    rownames(matrixTest) <- vapply(dxTestSVM[1], rownames,
                                FUN.VALUE = character(dim(dxTestSVM[[1]])[1]))
    return(matrixTest)
}

.svmTest <- function(svmMoldel, dxTest, rootOnto, variance) {
    dxTest <- as.data.frame(dxTest)
    matrixScore <- matrix(0, dim(dxTest)[1], 1)
    rownames(matrixScore) <- rownames(dxTest)
    nodeOnto <- svmMoldel$ontoTerm
    colnames(matrixScore) <- nodeOnto
    if (nodeOnto != rootOnto) {
        zeroValues <- svmMoldel$zeroValue
        if (length(zeroValues) != 0) dxTest <- dxTest[, -zeroValues]
        predOnto <- predict(svmMoldel, dxTest, decision.values = TRUE)
        matrixScore[, 1] <- 1 / (1 + exp(-2 * attr(predOnto, "decision.values")/
            variance[nodeOnto]))
    } else {
        matrixScore[, 1] <- 0.9999
    }
    matrixScore[which(matrixScore < 1.0e-06)] <- 1.0e-06
    matrixScore[which(matrixScore > 0.9999999)] <- 0.9999
    return(matrixScore)
}
