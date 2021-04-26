tableTPG <- function(att){
  x <- as.matrix(expand.grid(rep(list(c(0, 1)), att)))
  indexTable <- seq(dim(x)[1]/2+1, dim(x)[1]-1, 1)
  x <- x[order(x[, 1]), ]
  x <- cbind(x[-indexTable, ], array(1, dim(x)[1]-length(indexTable)))
  colnames(x) <- c("X_child",
                   sapply(seq_len(dim(x)[2]-2), FUN=function(att)
                     (paste("X_father", att, sep = "", collapse = ''))),
                   "TPG_rule")
  x <- x + 1
  return(x)
}
