\name{svmTrain}
\alias{svmTrain}
\title{
Binary SVM classification models for individual Ontology-term predictions
}
\description{
svmTrain delivers a set of binary SVM classifiers for different Ontology-terms.
}
\usage{
svmTrain(nodeGraph, tableOntoTerms, dxCharacterized, graphOnto,
        kernelSVM = "radial")
}

\arguments{
    \item{nodeGraph}{A character indicating a GO node label}

    \item{tableOntoTerms}{A binary matrix with `n' proteins (rows) by `m' Ontology node labels (columns).}

    \item{dxCharacterized}{A data frame with `n' protein coding genes (rows) by `f' features (columns).}

    \item{graphOnto}{A graphNEL graph with `m' Ontology node labels.}

    \item{kernelSVM}{The kernel used to calculate the variance (default: radial).}
}

\details{
Starting from sets of positively annotated protein sequences to different GO-terms in a GO subgraph, corresponding sets of negatively annotated protein sequences are computed using the inclusive separation policy proposed by Eisner et al. Training datasets for each GO-term are used to train binary Support Vector Machine (SVM) classifiers with a variety of kernel options.
}
\value{
    \item{svmTrain}{A list of objects of ``svm" class containing the fitted model.}
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
\code{\link[fgga]{svmOnto}}
}
\examples{
data(CfData)

mygraphGO <- as(CfData[["graphCfGO"]], "graphNEL")

# SVM training in four GO-terms
modelSVMs <- lapply(CfData[["nodesGO"]][1:4], FUN = svmTrain,
                    tableOntoTerms = CfData[["tableCfGO"]],
                    dxCharacterized = CfData[["dxCf"]],
                    graphOnto = mygraphGO, kernelSVM = "radial")

}

\keyword{svm}
