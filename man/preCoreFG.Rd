\name{preCoreFG}
\alias{preCoreFG}
\title{Transitive closure processing of a cross-ontology DAG}
\description{
preCoreFG ensures the transitive closure of inference paths -serial concatenation of relationships- in a cross-ontology DAG.}

\usage{
preCoreFG(ontoTerms, domains = "GO")
}

\arguments{
    \item{ontoTerms}{A vector with `m' cross-ontology node labels}

    \item{domains}{A string that indicates which subdomains or ontologies will be used. Values: ``GOBP", ``GOMF", ``GOCC", ``GOCC-PO", ``GOCC-ZFA", ``GOBP-HPO", ``GOMF-HPO", ``GOCC-HPO", ``GO-PO", ``GO-ZFA", ``GO-HPO", ``GO" (default, ``BP-MF-CC")}

}

\details{
Non-transitive relationships in cross-ontology DAG's may lead to non-transitive inference paths precluding the free propagation and consistency checking of ontology annotations. A transitive closure screening process over cross-ontology DAG's relationships is required before the construction of Forney Factor Graphs. Serial concatenation of relationships leading to non-transitive inference paths in a cross-ontology DAG are conformed by removing the most specific relationship.
}

\value{
A graphNEL graph with `m' node labels belong to ontologies used.
}

\references{
Spetale Flavio E., Arce D., Krsticevic F., Bulacio P. and Tapia E. “Consistent prediction of GO protein localization”. Scientific Report 7787(8), 2018
}
\author{
Flavio E. Spetale <spetale@cifasis-conicet.gov.ar>
}

\seealso{
\code{\link[fgga]{fgga2bipartite}}
}
\examples{

data(CfData)

myGOs <- c(CfData[["nodesGO"]], "GO:1902494", "GO:0032991", "GO:1990234",
            "GO:0005575")

\donttest{
# mygraphGO <- preCoreFG(myGOs, domains = "GOMF")
}}

\keyword{msgFFGA}
