\name{preCoreFG}
\alias{preCoreFG}
\title{Transitive closure processing of a GO DAG}
\description{
preCoreFG ensures the transitive closure of inference paths -serial concatenation of relationships- in a GO DAG.}

\usage{
preCoreFG(myGOnames, domains = "goPlus")
}

\arguments{
  \item{myGOnames}{A vector with `m' GO node labels}

  \item{domains}{A string that indicates which subdomains will be used. Values: ``BP-MF", ``MF-CC", ``BP-CC" or ``goPlus" (default, ``BP-MF-CC")}

}

\details{
Non-transitive relationships in GO DAG's may lead to non-transitive inference paths precluding the free propagation and consistency checking of GO annotations. A transitive closure screening process over GO DAG's relationships is required before the construction of Forney Factor Graphs. Serial concatenation of relationships leading to non-transitive inference paths in a GO DAG are conformed by removing the most specific relationship.
}
\value{
A graphNEL graph with `m' GO node labels.
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

\dontrun{
mygraphGO <- preCoreFG(myGOs, domains = "MF-CC")
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{msgFFGA}

