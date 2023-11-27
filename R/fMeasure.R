fMeasures <- function(target, predicted, cutoff = 0.5) {
    if ((nrow(target) != nrow(predicted)) || (ncol(target) != ncol(predicted)))
        stop ("Number of rows or columns do not match between target and
        predicted terms")
    predicted <- apply(predicted, c(1,2), function(x)
        {if (x >= cutoff) x <- 1 else x <- 0})
    ontoTerm <- colnames(target)
    measures <- c("Prec", "Recall", "Specif", "Fmeasure", "Acc", "nPositive")
    perfByTerms <- matrix(vapply(seq_len(ncol(target)),FUN = function(x,y,z)
        {.calculateMetrics(y[,x],z[,x])}, y=predicted, z=target,
        FUN.VALUE = numeric(6)), ncol(target), 6,
        dimnames = list(ontoTerm, measures), byrow = TRUE)

    avgTerm <- apply(perfByTerms, 2, mean)

    return(list(avgTerm=avgTerm, perfByTerms=perfByTerms))
}

fMeasuresByLevel <- function(target, predicted, graphOnto, cutoff = 0.5) {
    if (!missing(graphOnto)) {
        predicted <- apply(predicted, c(1,2), function(x)
        {if (x >= cutoff) x <- 1 else x <- 0})
        levels <- maxDistancegraphOnto(graphOnto)
        ontoTerm <- colnames(target)
        measures <- c("Prec", "Recall", "Specif", "Fmeasure", "Acc",
                    "nPositive")
        perfByTerms <- matrix(vapply(seq_len(ncol(target)),FUN = function(x,y,z)
            {.calculateMetrics(y[,x],z[,x])}, y=predicted, z=target,
            FUN.VALUE = numeric(6)), ncol(target), 6,
            dimnames = list(ontoTerm, measures), byrow = TRUE)

        gLevels <- array(sort(unique(levels)))
        perfByLevel <- matrix(0, length(gLevels), 6, dimnames=list(gLevels,
                                                                    measures))

        for (i in seq_len(length(gLevels))){
            if (gLevels[i]>0 && i==1) stop ("The graph does not have a root")
            if (gLevels[i]==0 && i==1) perfByLevel[i, ] <- apply(matrix(
                perfByTerms[which(levels == 0), ], nrow=1), 2, mean) else {
                indexLevels <- which(levels == gLevels[i])
                if (length(indexLevels) == 1) perfByLevel[i, ] <-
                    apply(t(perfByTerms[indexLevels, ]), 2, mean) else
                    perfByLevel[i, ] <- apply(
                        perfByTerms[indexLevels, ], 2, mean)
                }
            }
        } else stop ("Missing graphOnto")

    return(list(avgTerm=apply(perfByLevel, 2, mean), perfByLevel=perfByLevel))
}

.calculateMetrics <- function(termPred, termTarget) {

    if (length(termPred)!=length(termTarget))
        stop("The number of gene products per predicted term is not equal to the
            target")
    indexPosTerm <- which(termTarget > 0)
    indexNegTerm <- setdiff(seq_len(length(termTarget)), indexPosTerm)

    TP <- sum(termPred[indexPosTerm] > 0)
    FN <- sum(termPred[indexPosTerm] <= 0)
    TN <- sum(termPred[indexNegTerm] <= 0)
    FP <- sum(termPred[indexNegTerm] > 0)
    acc <- (TP+TN)/length(termTarget)

    if ((TP+FP) == 0) prec <- 0 else prec <- TP/(TP+FP)
    if ((TP+FN) == 0) recall <- 0 else recall <- TP/(TP+FN)
    if ((TN+FP) == 0)  specif <- 0 else specif <- TN/(TN+FP)
    if ((prec + recall) == 0) fMeasure <- 0 else
        fMeasure <- 2 *(prec * recall) / (prec + recall)

    measures <- c(prec,recall,specif,fMeasure,acc, length(indexPosTerm))
    names(measures) <- c("Prec", "Recall", "Specif", "Fmeasure", "Acc",
                        "nPositive")
    return (measures)
}

fHierarchicalMeasures <- function(target, predicted, graphOnto, cutoff = 0.5){
    if (dim(predicted)[2] != dim(target)[2])
        stop("The number of ontology term predicted is not equal to the target")

    commonId <- order(intersect(rownames(predicted), rownames(target)))
    target <- target[commonId, ]
    predicted <- predicted[commonId, ]
    samplesNoEval <- vector()
    gRoot <- leaves(graphOnto, "in")
    predicted <- apply(predicted, c(1,2), function(x)
        {if (x >= cutoff) x <- 1 else x <- 0})

    infoGraphBase <- list(graphOntoParents = list(),
                        graphOntoAncestors = list())
    ontoTerms <- nodes(graphOnto)
    infoGraphBase[["graphOntoParents"]] <- lapply(ontoTerms, FUN = function(x,y)
        {parents(x,y)}, y=graph_from_graphnel(graphOnto))
    names(infoGraphBase[["graphOntoParents"]]) <- ontoTerms
    infoGraphBase[["graphOntoAncestors"]] <- lapply(ontoTerms,
                    FUN = function(x,y) {ancestralSet(x,y)},
                    y=graph_from_graphnel(graphOnto))
    names(infoGraphBase[["graphOntoAncestors"]]) <- ontoTerms

    nSample <- j <- 0
    sumHR <- sumHP <- sumHF <- 0
    for (i in seq_len(dim(predicted)[1])){
        predTerms <- names(which(predicted[i,] == 1))
        targetTerms <- names(which(target[i,] == 1))
        if(length(predTerms) > 0 && length(targetTerms) > 2){
            nSample <- nSample+1
            HP <- .hPrecision(targetTerms, predTerms, graphOnto, gRoot,
                                infoGraphBase)
            HR <- .hRecall(targetTerms, predTerms, graphOnto, gRoot,
                                infoGraphBase)
            HF <- .hFmeasure(HP, HR)
            sumHP <- sumHP + HP
            sumHR <- sumHR + HR
            sumHF <- sumHF + HF
        } else {
            j<-j + 1
            samplesNoEval[j] <- rownames(predicted)[i]}
        }
    sumHP <- sumHP/nSample
    sumHR <- sumHR/nSample
    sumHF <- sumHF/nSample
    hMeasures <- list(sumHP, sumHR, sumHF, nSample, samplesNoEval)
    names(hMeasures) <- c("HP", "HR", "HF", "nSample", "noEvalSample")
    return(hMeasures)
}

.hPrecision<-function(target, predict, graphOnto, gRoot, graphOntoProperties){
    gTarget <- subGraph(target, graphOnto)
    gPred <- subGraph(predict, graphOnto)
    leafTarget <- leaves(gTarget, "out")
    leafPred <- leaves(gPred, "out")
    HP <- 0
    nLeafPred <- length(leafPred)
    nLeafTarget <- length(leafTarget)
    sumPred <- 0
    for (i in seq_len(nLeafPred)){
      parentsPred <- graphOntoProperties[["graphOntoAncestors"]][[
        leafPred[i]]]
      maxPred <- array(0,nLeafTarget)
      for (j in seq_len(nLeafTarget)){
        if(leafTarget[j] != gRoot){
          parentsTarget <- graphOntoProperties[["graphOntoAncestors"]][[
                                              leafTarget[j]]]
        }
        else parentsTarget <- gRoot
        maxPred[j] <- length(intersect(parentsPred, parentsTarget)) /
          length(parentsPred)
      }
      sumPred <- max(maxPred) + sumPred
    }
    HP <- sumPred / nLeafPred
    return(HP)
}

.hRecall<-function(target, predict, graphOnto, gRoot, graphOntoProperties){
    gTarget <- subGraph(target, graphOnto)
    gPred <- subGraph(predict, graphOnto)
    leafTarget <- leaves(gTarget, "out")
    leafPred <- leaves(gPred, "out")
    HR <- 0
    nLeafPred <- length(leafPred)
    nLeafTarget <- length(leafTarget)
    sumPred <- 0
    for (i in seq_len(nLeafTarget)){
        parentsTarget <- graphOntoProperties[["graphOntoAncestors"]][[
        leafTarget[i]]]
        maxPred <- array(0,nLeafPred)
        for (j in seq_len(nLeafPred)){
            if(leafPred[j] != gRoot){
                parentsPred <- graphOntoProperties[["graphOntoAncestors"]][[
                leafPred[j]]]
                }
        else parentsPred <- gRoot
            maxPred[j] <- length(intersect(parentsPred, parentsTarget)) /
            length(parentsTarget)
        }
        sumPred <- max(maxPred) + sumPred
    }
    HR <- sumPred / nLeafTarget
    return(HR)
}

.hFmeasure<-function(HP, HR){
    if (HP==0 && HR == 0) return (0) else return(2*HP*HR / (HP + HR))
}
