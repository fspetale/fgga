\name{svmOnto}
\alias{svmOnto}
\title{
Ontology-term predictions by binary SVM classifiers
}
\description{
svmOnto delivers soft Ontology-term predictions based on binary SVM classification models.
}
\usage{
svmOnto(svmMoldel, dxCharacterized, rootNode, varianceSVM)
}
\arguments{

    \item{svmMoldel}{A list of object of class ``svm" created by svm.}

    \item{dxCharacterized}{A data frame with `n' protein coding genes (rows) by `f' features (columns).}

    \item{rootNode}{A character indicating the root of the graph.}

    \item{varianceSVM}{A vector named with the variance of cross-Ontology node labels.}
}
\details{
Binary SVM predictions are supplemented with their corresponding margins. These margins are used to model the additive zero-mean Gaussian learning noise that corrupts ideal but hidden Ontology-term predictions. These ideal predictions are embedded in hidden variable nodes of the Forney Factor Graph.
}
\value{
    \item{svmOnto}{A named vector of predicted values for a protein sequence.}
}
\references{
Chang, Chih-Chung and Lin, Chih-Jen: LIBSVM: a library for Support Vector Machines
\url{http://www.csie.ntu.edu.tw/~cjlin/libsvm}

Eisner R, Poulin B, Szafron D, Lu P, Greiner R. Improving protein function prediction using the hierarchical structure of the Gene Ontology. In: Proc. IEEE CIBCB; 2005. p. 1–1

Spetale FE, Tapia E, Krsticevic F, Roda F, Bulacio P (2016). A Factor Graph Approach to Automated GO Annotation. PLOS ONE 11(1): e0146986
}
\author{
Flavio E. Spetale, Pilar Bulacio and Javier Murillo <spetale@cifasis-conicet.gov.ar>
}

\seealso{
\code{\link[fgga]{svmTrain}}
}
\examples{
data(CfData)

mygraphGO <- as(CfData[["graphCfGO"]], "graphNEL")

modelSVMs <- lapply(CfData[["nodesGO"]][1:4], FUN = svmTrain,
                    tableOntoTerms = CfData[["tableCfGO"]],
                    dxCharacterized = CfData[["dxCf"]],
                    graphOnto = mygraphGO, kernelSVM = "radial")

rootGO <- leaves(mygraphGO, "in")

varianceGOs <- CfData[["varianceGOs"]]

# SVM testing in four GO-terms
dxTestCharacterized <- CfData[["dxCf"]][
                                sample(1:dim(CfData[["dxCf"]])[1], 20), ]

matrixGOTest <- svmOnto(svmMoldel = modelSVMs,
                        dxCharacterized = dxTestCharacterized,
                        rootNode = rootGO, varianceSVM = varianceGOs)
}

\keyword{svm}
