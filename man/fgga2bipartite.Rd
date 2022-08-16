\name{fgga2bipartite}
\alias{fgga2bipartite}
\title{
Forney Factor Graph model
}
\description{
fgga2bipartite builds a Forney Factor Graph from a FGGA model.
}
\usage{
fgga2bipartite(graphOnto)
}

\arguments{
    \item{graphOnto}{A graphNEL graph with `m' cross-Ontology node labels.}
}

\details{
The \bold{Gene Ontology} (GO) is structured as a directed acyclic graph (DAG) with nodes (GO-terms) representing gene functions and edges characterizing relationships between nodes. A variety of relationships are possible (currently 8). To compute GO-term predictions perfectly aware of GO-term relationships, a  Forney Factor Graph is required. Hence, GO-terms are mapped to binary variable nodes, and relationships to logical factor nodes.
}

\value{
A binary matrix with \emph{2*m} rows by \emph{2*m-1} columns where \emph{m} is the quantity of cross-Ontology node labels.
}

\references{
F. Spetale, J. Murillo, E. Tapia, D. Arce, S. Ponce, and P. Bulacio, “Formal modeling of gene ontology annotation predictions based on factor graphs,” Journal of Physics: Conference Series, vol. 705, no. 1, p. 012001, 2016.

Spetale F.E., Tapia E., Krsticevic F., Roda F. and Bulacio P. “A Factor Graph Approach to Automated GO Annotation”. PLoS ONE 11(1): e0146986, 2016.

Spetale Flavio E., Arce D., Krsticevic F., Bulacio P. and Tapia E. “Consistent prediction of GO protein localization”. Scientific Report 7787(8), 2018
}

\author{
Flavio E. Spetale <spetale@cifasis-conicet.gov.ar>
}

\examples{
data(CfData)

graphGO <- as(CfData$graphCfGO, "graphNEL")
fgga2bipartite(graphGO)
}
