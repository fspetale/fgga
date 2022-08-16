library(fgga)
test_graph <- function(){
  V <- LETTERS[1:10]
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

test_msg <- function() {

  g <- test_graph()
  mF <- fgga2bipartite(g)
  mT <- matrix(runif(20, 0, 1), 2, 10)
  colnames(mT) <- nodes(g)
  rownames(mT) <- sapply(1:2, FUN = function(x) paste0("IDT", x, sep = ""))
  mT[ , 1] <- 0.9999

  mR <- t(apply(mT, MARGIN = 1, FUN = msgFGGA, matrixFGGA = mF,
                            graphOnto = g, tmax = 50, epsilon = 0.001))

  checkEquals(length(mR), 20)
  checkIdentical(rownames(mR)[1], "IDT1", "Incorrect Test")
  checkEqualsNumeric(sum(mR[,1])/2, 1, tolerance = 0.01)
}
