\name{fgga}
\alias{fgga}

\title{
Factor Graph GO Annotation model
}
\description{
A hierarchical graph-based machine learning model for the consistent GO
annotation of protein coding genes.
}
\usage{
fgga(graphGO, tableGOs, dxCharacterized, dxTestCharacterized,
    kFold, kernelSVM, tmax, epsilon)
}

\arguments{
    \item{graphGO}{A graphNEL graph with `m' GO node labels.}

    \item{tableGOs}{A binary matrix with `n' proteins (rows) by `m' GO node labels (columns).}

    \item{dxCharacterized}{A data frame with `n' proteins (rows) by `f' features (columns).}

    \item{dxTestCharacterized}{A data frame with `k' proteins (rows) by `f' features (columns).}

    \item{kFold}{An integer for the number of folds.}

    \item{kernelSVM}{The kernel used to calculate the variance (default: radial).}

    \item{tmax}{An integer indicating the maximum number of iterations (default: 200).}

    \item{epsilon}{An integer that represents the convergence criteria (default: 0.001).}
}

\details{
The \bold{FGGA model} is built in two main steps. In the first step, a core Factor Graph (FG) modeling hidden GO-term predictions and relationships is created. In the second step, the FG is enriched with nodes modeling observable GO-term predictions issued by \bold{binary SVM classifiers}. In addition, probabilistic constraints modeling learning gaps between hidden and observable GO-term predictions are introduced. These gaps are assumed to be independent among GO-terms, locally additive with respect to observed predictions, and zero-mean Gaussian. \bold{FGGA predictions} are issued by the native iterative \bold{message passing algorithm} of factor graphs.
}

\value{
A named matrix with `k' protein coding genes (rows) by `m' GO node labels (columns) where each element indicates a probabilistic prediction value.
}

\references{
Spetale F.E., Tapia E., Krsticevic F., Roda F. and Bulacio P. “A Factor Graph Approach to Automated GO Annotation”. PLoS ONE 11(1): e0146986, 2016.

Spetale Flavio E., Arce D., Krsticevic F., Bulacio P. and Tapia E. “Consistent prediction of GO protein localization”. Scientific Report 7787(8), 2018
}

\author{
Flavio E. Spetale and Elizabeth Tapia <spetale@cifasis-conicet.gov.ar>
}

\seealso{
\code{\link[fgga]{fgga2bipartite}}, \code{\link[fgga]{sumProduct}}, \code{\link[fgga]{svmGO}}
}
\examples{
data(CfData)

mygraphGO <- as(CfData[["graphCfGO"]], "graphNEL")

dxCfTestCharacterized <- CfData[["dxCf"]][CfData[["indexGO"]]$indexTest[1:2], ]

myTableGO <- CfData[["tableCfGO"]][
                    CfData[["indexGO"]]$indexTrain[1:300], ]

dataTrain <- CfData[["dxCf"]][
                    CfData[["indexGO"]]$indexTrain[1:300], ]

fggaResults <- fgga(graphGO = mygraphGO,
                tableGOs = myTableGO, dxCharacterized = dataTrain,
                dxTestCharacterized = dxCfTestCharacterized, kFold = 2,
                tmax = 50, epsilon = 0.05)

}

\keyword{fgga2bipartite}
\keyword{msgFGGA }
