maxDistancegraphOnto <- function(graphOnto){
    terms <- nodes(graphOnto)
    matrixDistgraphOnto <- matrix(0,length(terms),2)
    matrixDistgraphOnto[1, 1] <- gRoot <- leaves(graphOnto, "in")
    matrixDistgraphOnto <- matrix("list", length(terms))
    termUsed <- array(0, length(terms))
    names(matrixDistgraphOnto) <- names(termUsed) <- terms
    gLeaves <- leaves(graphOnto, "out")
    pNodes <- terms
    while(length(pNodes)>0){
        pathParents <- parents(pNodes[1], graphOnto)
        if (length(pNodes) == length(terms)){
            pNodes <- setdiff(terms, gRoot)
            matrixDistgraphOnto[gRoot] <- gRoot
            termUsed[gRoot] <- 1
        }
        if (length(pathParents) > 0){
            allNodes <- intersect(pathParents, names(
                termUsed[which(termUsed > 0)]))
            if (length(allNodes) == length(pathParents)){
                matrixDistgraphOnto[pNodes[1]] <- as.data.frame(c(pNodes[1],
                as.character(unlist(matrixDistgraphOnto[pathParents[
                which.max(termUsed[pathParents])]]))))
            termUsed[pNodes[1]] <- termUsed[pathParents[
                which.max(termUsed[pathParents])]] + length(pNodes[1])
            pNodes <- pNodes[-1]
            } else {
            nodesUnused <- setdiff(pathParents,names(termUsed[
                which(termUsed > 0)]))
            pLastNode <- which(pNodes == nodesUnused[length(nodesUnused)])
            if (pLastNode == length(pNodes)){
                pNodes <- c(pNodes[seq(2, pLastNode)], pNodes[1])
            } else pNodes<-c(pNodes[seq(2, pLastNode)], pNodes[1],
                            pNodes[seq((pLastNode + 1), length(pNodes))])
            }
        }
    }
    return(termUsed=termUsed-1)
}
