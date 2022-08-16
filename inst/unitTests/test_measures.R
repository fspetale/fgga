library(fgga)
test_graph <- function(){
    V <- c("GO:0000166", "GO:0003674", "GO:0003676", "GO:0003677", "GO:0003824",
          "GO:0004888", "GO:0005102", "GO:1901363", "GO:1902494", "GO:1990234")
    edL <- vector("list", length=length(V))
    names(edL) <- V
    edL[[1]] <- list(edges=c(2, 3, 7, 5))
    edL[[2]] <- list(edges=c(4, 6))
    edL[[3]] <- list(edges=c(6, 8))
    edL[[4]] <- list(edges=c(6, 9))
    edL[[5]] <- list(edges=c(8))
    edL[[6]] <- list(edges=c(8, 9, 10))
    edL[[7]] <- list(edges=c(8, 5))
    edL[[8]] <- list(edges=c())
    edL[[9]] <- list(edges=c())
    g <- graph::graphNEL(nodes=V, edgeL=edL, edgemode="directed")
    return(g)
}

test_maxDistGraph <- function(){
    g <- test_graph()
    distGraph <- maxDistancegraphOnto(g)
    checkEquals(length(distGraph), length(nodes(g)))
    checkTrue(is.numeric(distGraph[1]))
}

test_measure <- function(){
    g <- test_graph()
    dxReal <- matrix(rbinom(360, 1, 0.4), 36, 10)
    colnames(dxReal) <- nodes(g)
    dxReal[, 1] <- 1
    dxTest <- matrix(runif(360, 0, 1), 36, 10)
    colnames(dxTest) <- nodes(g)
    dxTest[, 1] <- 1
    xTest <- fMeasures(dxReal, dxTest, cutoff = 0.5)
    checkIdentical(xTest[["perfByTerms"]][1, 6], 36, "Incorrect Test")
    checkEqualsNumeric(xTest[["perfByTerms"]][1, 1], 1, tolerance = 0.01)

    rownames(dxTest) <- rownames(dxReal) <- paste0("p", seq_len(36), sep="")
    hTest <- fHierarchicalMeasures(dxReal, dxTest, g, cutoff = 0.5)
    checkEquals(length(hTest), 5)
    checkTrue(is.numeric(hTest[["HP"]]))
}
