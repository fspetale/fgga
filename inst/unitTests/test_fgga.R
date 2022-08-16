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

test_fgga <- function(){
    g <- test_graph()
    g <- as(g, "matrix")
    g <- rbind(g, array(0,10))
    g <- cbind(g, array(0,11))
    colnames(g)[11] <- rownames(g)[11] <- "Ontology:FES"
    g["Ontology:FES", "GO:0000166"] <- 1
    g <- as(g, "graphNEL")
    dx1 <- matrix(runif(1800, 0, 1), 360, 5)
    dx2 <- matrix(rbinom(3960, 1, 0.4), 360, 11)
    dx2[, 11] <- rep(1, 360)
    dx3 <- matrix(runif(10, 0, 1), 2, 5)
    colnames(dx1) <- colnames(dx3) <- c("AA", "BB", "CC", "DD", "EE")
    rownames(dx1) <- rownames(dx2) <- sapply(seq_len(360),
                    FUN = function(x) paste0("ID", x, sep = ""))
    rownames(dx3) <- sapply(seq_len(2),
                    FUN = function(x) paste0("IDT", x, sep = ""))
    colnames(dx2) <- nodes(g)

    mR <- fgga(graphOnto = g, tableOntoTerms = dx2, dxCharacterized = dx1,
                        dxTestCharacterized = dx3, kFold = 2,
                        tmax = 50, epsilon = 0.05)

    checkEquals(length(mR), 22)
    checkIdentical(rownames(mR)[2], "IDT2", "Incorrect Test")
    checkEqualsNumeric(sum(mR[,11])/2, 1, tolerance = 0.01)

}
