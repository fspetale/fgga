\name{varianceOnto}
\alias{varianceOnto}
\title{
The variance of the gaussian learning noise at individual Ontology-terms
}
\description{
varianceOnto estimates the variance of gaussian distributions modeling the additive learning noise that corrupts ideal Ontology-term predictions.
}
\usage{
varianceOnto(tableOntoTerms, dxCharacterized, kFold, graphOnto, rootNode,
            kernelSVM = "radial")
}
\arguments{
    \item{tableOntoTerms}{A binary matrix with `n' protein coding genes (rows) by `m' cross-Ontology node labels (columns).}

    \item{dxCharacterized}{A data frame with `n' protein coding genes (rows) by `f' features (columns).}

    \item{kFold}{An integer for the number of folds.}

    \item{graphOnto}{A graphNEL graph with `m' cross-Ontology node labels.}

    \item{rootNode}{A character indicating the root of the graph.}

    \item{kernelSVM}{The kernel used to calculate the variance (default: radial).}
}
\details{
Under the assumption of symmetrical (Gaussian) conditional probability distributions for observable variable node predictions \eqn{y_i} over a hidden variable node annotations \eqn{x_i}, variances \eqn{\eta_i} can be estimated using a validation dataset of positively annotated samples. \if{latex}{\out{Let \emph{D} be a validation dataset with $L^{+}$ positively annotated samples
\deqn{\hat{\eta_i} = 1/(L^{+} - 1) * \sum_{i=1}^L (x_i - y_i)}

where $x_i = 1$ is the positive annotation of the l-th data sample to the ith Ontology-term and $y_i$  is the corresponding real-valued classifier (SVM) prediction.}}
}

\value{
A vector named with the variance of each cross-Ontology node.}

\references{
Spetale FE, Tapia E, Krsticevic F, Roda F, Bulacio P (2016). A Factor Graph Approach to Automated GO Annotation. PLOS ONE 11(1): e0146986
}
\author{
Flavio E. Spetale <spetale@cifasis-conicet.gov.ar>
}

\examples{
data(CfData)

mygraphGO <- as(CfData[["graphCfGO"]], "graphNEL")

rootGO <- leaves(mygraphGO, "in")

mygraphGO <- subGraph(c("GO:0140110", "GO:0098772", "GO:0003674"), mygraphGO)

myTableGO <- CfData[["tableCfGO"]][
                    CfData[["indexGO"]]$indexTrain,
                    c("GO:0140110", "GO:0098772", "GO:0003674")]

varianceGOs <- varianceOnto(tableOntoTerms = myTableGO,
                            dxCharacterized = CfData[["dxCf"]],
                            kFold = 2, graphOnto = mygraphGO,
                            rootNode = rootGO, kernelSVM = "radial")
}

\keyword{variance}
