\name{fgga-package}
\alias{fgga-package}
\docType{package}
\title{FGGA: Factor Graph Gene ontology Annotation.}
\description{
FGGA is a graph-based machine learning approach for the automated and consistent GO annotation of protein coding genes. The input is a set of GO-term annotated protein coding genes previously characterized in terms of a fixed number of user-defined features, including the presence/absence of PFAM domains, physical-chemical properties, presence of signal peptides, among others. The set of GO-terms defines the output GO subgraph. A hierarchical ensemble (SVMs) machine learning model is generated. This model can be used to predict the GO subgraph annotations of uncharacterized protein coding genes. Individual GO-term annotations are accompanied by maximum a posteriori probability estimates issued by the native message passing algorithm of factor graphs.
}
\details{
\tabular{ll}{
Package: \tab fgga\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2021-04-06\cr
License: \tab GPL 3.0\cr
LazyLoad: \tab yes\cr
}
}
\author{
\emph{Flavio E. Spetale}, \emph{Javier Murillo} and \emph{Elizabeth Tapia}

BioInformatics

Cifasis-Conicet

\email{spetale@cifasis-conicet.gov.ar}

Maintainer:
\emph{Flavio E. Spetale}
}

\references{
Spetale F.E., et al. \bold{A Factor Graph Approach to Automated GO Annotation}. \emph{PLoS ONE}  (2016). \url{ https://doi.org/10.1371/journal.pone.0146986}.

Spetale Flavio E., et al. \bold{Consistent prediction of GO protein localization}. \emph{Scientific Report} (2018). \url{https://doi.org/10.1038/s41598-018-26041-z}.

}
\keyword{ package }
\seealso{
\code{\link[fgga]{fgga}}, \code{\link[fgga]{fgga2bipartite}}, \code{\link[fgga]{sumProduct}}, \code{\link[e1071]{svm}}
}
\examples{
# Loading example dataset
data(CfData)

# Selecting training and testing data
indexTrain <- sample(1:dim(CfData[["dxCf"]])[1], 4500)

dxCfTrainCharacterized <- CfData[["dxCf"]][indexTrain, ]

indexTest <- setdiff(seq_len(dim(CfData[["dxCf"]])[1]), indexTrain)

dxCfTestCharacterized <- CfData[["dxCf"]][indexTest, ]

\donttest{
# Applying the FGGA algorithm
fggaResults <- fgga(graphGO = CfData[["graphCfGO"]],
                    tableGOs = CfData[["tableCfGO"]],
                    dxCharacterized = dxCfTrainCharacterized,
                    dxTestCharacterized = dxCfTestCharacterized,
                    kFold = 5, tmax = 50, epsilon = 0.001)

}}
