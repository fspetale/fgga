\name{sumProduct}
\alias{msgFGGA}
\title{
Message passing algorithm between nodes of the Forney Factor Graph
}
\description{
msgFGGA operates in Forney Factor Graphs and computes approximate maximum a posteriori (MAP) estimates of hidden Ontology variable nodes (Ontology-terms).
}

\usage{
msgFGGA(matrixFGGA, obsValueOntoTerms, graphOnto, tmax = 200,
        epsilon = 0.001)
}

\arguments{
    \item{matrixFGGA}{A binary matrix with FGGA model of the class `fgga.'}

    \item{obsValueOntoTerms}{A named vector with `m' probabilistic prediction values for a protein coding gene.}

    \item{graphOnto}{A graphNEL graph with `m' Ontology node labels.}

    \item{tmax}{An integer indicating the maximum number of iterations (defualt: 200).}

    \item{epsilon}{An integer that represents the convergence criteria (default: 0.001)}
}

\details{
Starting from Ontology-term predictions at observable variable nodes, probability distribution functions modelling the learning noise of individual Ontology-terms, a user-defined number of iterations (maximum 200), a user-defined threshold for the convergence of predictions (maximum 0.001), and the structure of the Forney Factor Graph, the \bold{msgFGGA} delivers approximate maximum a posteriori (MAP) estimates of hidden GO variable nodes (GO-terms).
}

\value{
A named vector with `m' consistent probabilistic predictions for a protein coding genes.
}
\references{
Kschischang FR, Frey BJ, Loeliger H.-A. Factor graphs and the sum-product algorithm. IEEE Trans. Inf. Theor. 47, 498–519 (2001).

Yedidia JS. Message-passing algorithms for inference and optimization. Journal of Statistical Physics 145, 860–890 (2011).

Spetale FE, Tapia E, Krsticevic F, Roda F, Bulacio P (2016). A Factor Graph Approach to Automated GO Annotation. PLOS ONE 11(1): e0146986
}
\author{
Flavio E. Spetale and Elizabeth Tapia <spetale@cifasis-conicet.gov.ar>
}

\seealso{
\code{\link[fgga]{tableTPG}}
}
\examples{
data(CfData)

mygraphGO <- as(CfData[["graphCfGO"]], "graphNEL")

myTableGO <- CfData[["tableCfGO"]][
                    CfData[["indexGO"]]$indexTrain[1:500], ]

modelSVMs <- lapply(CfData[["nodesGO"]], FUN = svmTrain,
                    tableOntoTerms = myTableGO,
                    dxCharacterized = CfData[["dxCf"]],
                    graphOnto = mygraphGO, kernelSVM = "radial")

rootGO <- leaves(mygraphGO, "in")

varianceGOs <- CfData[["varianceGOs"]]

dxTestCharacterized <- CfData[["dxCf"]][
                                sample(1:dim(CfData[["dxCf"]])[2], 2), ]

matrixGOTest <- svmOnto(svmMoldel = modelSVMs,
                        dxCharacterized = dxTestCharacterized,
                        rootNode = rootGO, varianceSVM = varianceGOs)

modelFGGA <- fgga2bipartite(mygraphGO)

matrixFGGATest <- t(apply(matrixGOTest, MARGIN = 1, FUN = msgFGGA,
                    matrixFGGA = modelFGGA, graphOnto= mygraphGO,
                    tmax = 50, epsilon = 0.1))

}

\keyword{msgFGGA}
