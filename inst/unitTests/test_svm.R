library(fgga)
test_graph <- function(){
    V <- sapply(seq_len(9), FUN = function(x) {
        if (x == 6) "Ontology:FES" else paste0("GO:000", x, sep = "")})
    g <- matrix(0, 9, 9)
    colnames(g) <- rownames(g) <- V
    g[1, 9] <- 1
    g[2, 9] <- 1
    g[3, c(2, 4, 7, 8)] <- 1
    g[4, 1] <- 1
    g[5, 1] <- 1
    g[6, 3] <- 1
    g[7, 5] <- 1
    g[8, 5] <- 1
    g <- as(g, "graphNEL")
    return (g)
}

test_svm <- function() {
    dx1 <- matrix(runif(360, 0, 1), 40, 9)
    dx3 <- matrix(runif(18, 0, 1), 2, 9)
    dx2 <- matrix(rbinom(360, 1, 0.3), 40, 9)
    dx2[, 6] <- rep(1, 40)
    colnames(dx1) <- colnames(dx2) <- colnames(dx3) <- V <- sapply(seq_len(9),
        FUN = function(x)
            { if (x == 6) "Ontology:FES" else paste0("GO:000", x, sep = "")})
    rownames(dx1) <- rownames(dx2) <-
        sapply(seq_len(40), FUN = function(x) paste0("ID", x, sep = ""))
    rownames(dx3) <- sapply(seq_len(2),
                            FUN = function(x) paste0("IDT", x, sep = ""))

    g <- test_graph()

    svmM <- lapply(V, FUN = svmTrain, tableOntoTerms = dx2, dxCharacterized = dx1,
                    graphOnto = g, kernelSVM = "radial")

    root <- leaves(g, "in")

    vx <- array(runif(9, 0, 1), 9)
    names(vx) <- colnames(dx1)
    classSVM <- vapply(svmM, class, FUN.VALUE = character(1))
    sinSVM <- which(classSVM=="list")
    if (length(sinSVM)>1){
        sinSVM <- setdiff(sinSVM, 6)
        svmM <- svmM[-sinSVM]
    }

    mT <- svmOnto(svmM, dx3, root, vx)

    checkTrue(is.numeric(mT))
    checkEquals(mT[1, 6], 0.9999)
}
