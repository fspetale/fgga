fgga2bipartite <- function(graphOnto) {
    matrixOnto <- as(graphOnto, "matrix")
    matrixOnto <- matrixOnto[, order(colnames(matrixOnto))]
    matrixOnto <- matrixOnto[order(rownames(matrixOnto)), ]
    nOnto <- dim(matrixOnto)[1]
    matrixFGGA <- matrix(0, nOnto * 2, nOnto * 2 - 1)
    rownames(matrixFGGA) <- c(rownames(matrixOnto), vapply(rownames(matrixOnto),
        FUN = function(att) {(paste("Obs_", att, sep = "", collapse = ""))},
        FUN.VALUE = character(1)))
    k <- vapply(seq_len(dim(matrixOnto)[1] - 1), FUN = function(att) {paste(
        "f", att, sep = "", collapse = "")}, FUN.VALUE = character(1))
    colnames(matrixFGGA) <- c(k, vapply(seq_len(dim(matrixOnto)[1]),
        FUN = function(att) {paste("g", att, sep = "", collapse = "")},
        FUN.VALUE = character(1)))
    linksFunctions <- apply(matrixOnto, MARGIN = 2,
                            FUN = function(x) (which(x == 1)))
    k <- 1
    for (i in seq_len(nOnto)) {
        matrixFGGA[c(i, i + nOnto), i + nOnto - 1] <- 1
        if (length(linksFunctions[[i]]) != 0) {
            matrixFGGA[c(i, linksFunctions[[i]]), k] <- 1
            k <- k + 1
        }
    }
    class(matrixFGGA) <- "fgga"
    return(matrixFGGA)
}
