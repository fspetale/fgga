\name{fgga}
\alias{fgga}

\title{
Factor Graph Cross-Ontology Annotation model
}
\description{
A hierarchical graph-based machine learning model for the consistent GO, PO, ZFA, HPO
annotation of protein coding genes.
}
\usage{
fgga(graphOnto, tableOntoTerms, dxCharacterized, dxTestCharacterized,
    kFold, kernelSVM, tmax, epsilon)
}

\arguments{
    \item{graphOnto}{A graphNEL graph with `m' Ontology node labels.}

    \item{tableOntoTerms}{A binary matrix with `n' proteins (rows) by `m' Ontology node labels (columns).}

    \item{dxCharacterized}{A data frame with `n' proteins (rows) by `f' features (columns).}

    \item{dxTestCharacterized}{A data frame with `k' proteins (rows) by `f' features (columns).}

    \item{kFold}{An integer for the number of folds.}

    \item{kernelSVM}{The kernel used to calculate the variance (default: radial).}

    \item{tmax}{An integer indicating the maximum number of iterations (default: 200).}

    \item{epsilon}{A real value less than 1 that represents the convergence criteria (default: 0.001).}
}

\details{
The \bold{FGGA model} is built in two main steps. In the first step, a core Factor Graph (FG) modeling hidden Ontology-term predictions and relationships is created. In the second step, the FG is enriched with nodes modeling observable Ontology-term predictions issued by \bold{binary SVM classifiers}. In addition, probabilistic constraints modeling learning gaps between hidden and observable Ontology-term predictions are introduced. These gaps are assumed to be independent among Ontology-terms, locally additive with respect to observed predictions, and zero-mean Gaussian. \bold{FGGA predictions} are issued by the native iterative \bold{message passing algorithm} of factor graphs.
}

\value{
A named matrix with `k' protein coding genes (rows) by `m' cross-Ontology node labels (columns) where each element indicates a probabilistic prediction value.
}

\references{
Spetale F.E., Tapia E., Krsticevic F., Roda F. and Bulacio P. “A Factor Graph Approach to Automated GO Annotation”. PLoS ONE 11(1): e0146986, 2016.

Spetale Flavio E., Arce D., Krsticevic F., Bulacio P. and Tapia E. “Consistent prediction of GO protein localization”. Scientific Report 7787(8), 2018
}

\author{
Flavio E. Spetale and Elizabeth Tapia <spetale@cifasis-conicet.gov.ar>
}

\seealso{
\code{\link[fgga]{fgga2bipartite}}, \code{\link[fgga]{sumProduct}}, \code{\link[fgga]{svmOnto}}
}
\examples{
data(CfData)

mygraphGO <- as(CfData[["graphCfGO"]], "graphNEL")

dxCfTestCharacterized <- CfData[["dxCf"]][CfData[["indexGO"]]$indexTest[1:2], ]

myTableGO <- CfData[["tableCfGO"]][
                    CfData[["indexGO"]]$indexTrain[1:300], ]

dataTrain <- CfData[["dxCf"]][
                    CfData[["indexGO"]]$indexTrain[1:300], ]

fggaResults <- fgga(graphOnto = mygraphGO,
                tableOntoTerms = myTableGO, dxCharacterized = dataTrain,
                dxTestCharacterized = dxCfTestCharacterized, kFold = 2,
                tmax = 50, epsilon = 0.05)

}

\keyword{fgga2bipartite}
\keyword{msgFGGA}
