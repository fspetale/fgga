\name{fgga-package}
\alias{fgga-package}
\docType{package}
\title{FGGA: Factor Graph Gene ontology Annotation.}
\description{
FGGA is a graph-based machine learning approach for the automated and consistent GO, PO, HPO and ZFA annotation of protein coding genes. The input is a set of ontoligical-terms annotated protein coding genes previously characterized in terms of a fixed number of user-defined features, including the presence/absence of PFAM domains, physical-chemical properties, presence of signal peptides, among others. The set of ontoligical terms defines the output cross-ontology subgraph. A hierarchical ensemble (SVMs) machine learning model is generated. This model can be used to predict the cross-ontology subgraph annotations of uncharacterized protein coding genes. Individual ontoligical-term annotations are accompanied by maximum a posteriori probability estimates issued by the native message passing algorithm of factor graphs.
}

\author{
\emph{Flavio E. Spetale}, \emph{Paolo Cacchiarrelli} and \emph{Elizabeth Tapia}

BioInformatics

Cifasis-Conicet

\email{spetale@cifasis-conicet.gov.ar}

Maintainer:
\emph{Flavio E. Spetale}
}

\references{
Spetale F.E., et al. \bold{A Factor Graph Approach to Automated GO Annotation}. \emph{PLoS ONE} (2016). \url{ https://doi.org/10.1371/journal.pone.0146986}.

Spetale Flavio E., et al. \bold{Consistent prediction of GO protein localization}. \emph{Scientific Report} (2018). \url{https://doi.org/10.1038/s41598-018-26041-z}.

}
\keyword{ package }
\seealso{
\code{\link[fgga]{fgga}}, \code{\link[fgga]{fgga2bipartite}}, \code{\link[fgga]{sumProduct}}, \code{\link[e1071]{svm}}
}
