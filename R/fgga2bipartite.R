fgga2bipartite <- function(graphGO) {
    matrixGO <- as(graphGO, "matrix")
    matrixGO <- matrixGO[, order(colnames(matrixGO))]
    matrixGO <- matrixGO[order(rownames(matrixGO)), ]
    nGO <- dim(matrixGO)[1]
    matrixFGGA <- matrix(0, nGO * 2, nGO * 2 - 1)
    rownames(matrixFGGA) <- c(rownames(matrixGO), vapply(rownames(matrixGO),
        FUN = function(att) {(paste("Obs_", att, sep = "", collapse = ""))},
        FUN.VALUE = character(1)))
    k <- vapply(seq_len(dim(matrixGO)[1] - 1), FUN = function(att) {paste(
        "f", att, sep = "", collapse = "")}, FUN.VALUE = character(1))
    colnames(matrixFGGA) <- c(k, vapply(seq_len(dim(matrixGO)[1]),
        FUN = function(att) {paste("g", att, sep = "", collapse = "")},
        FUN.VALUE = character(1)))
    linksFunctions <- apply(matrixGO, MARGIN = 2,
                            FUN = function(x) (which(x == 1)))
    k <- 1
    for (i in seq_len(nGO)) {
        matrixFGGA[c(i, i + nGO), i + nGO - 1] <- 1
        if (length(linksFunctions[[i]]) != 0) {
            matrixFGGA[c(i, linksFunctions[[i]]), k] <- 1
            k <- k + 1
        }
    }
    class(matrixFGGA) <- "fgga"
    return(matrixFGGA)
}
