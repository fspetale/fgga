\name{maxDistancegraphOnto}
\alias{maxDistancegraphOnto}

\title{
Maximum distance for a graph
}
\description{
Computes the maximum distance from any node to the root of the graph
}
\usage{
maxDistancegraphOnto(graphOnto)
}

\arguments{
    \item{graphOnto}{A graphNEL graph with `m' Ontology node labels.}
}

\details{
This function computes a distance matrix for a graph
}

\value{
Named numeric array containing the distance from any node to the root.
}

\seealso{
\code{\link[fgga]{fMeasure}}
}
\author{
Flavio E. Spetale <spetale@cifasis-conicet.gov.ar>
}

\examples{
data(CfData)

mygraphGO <- as(CfData[["graphCfGO"]], "graphNEL")

maxDistancegraphOnto(mygraphGO)
}
