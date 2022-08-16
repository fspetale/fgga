\name{tableTPG}
\alias{tableTPG}

\title{
Valid configurations for hidden variable nodes in a Forney Factor Graph
}
\description{
tableTPG provides valid configurations of hidden variable nodes at logical function nodes in a Forney Factor Graph under the True Path Graph (TPG) constraint.
}
\usage{
tableTPG(att)
}

\arguments{
    \item{att}{An integer indicating the number of cross-Ontology nodes involved}
}
\details{
Valid configurations of hidden variable nodes at logical function nodes enable messaging passing across the Forney Factor Graph. The TPG constraint is defined as: “If the child Ontology node describes the protein, then all its parent terms must also apply to that protein; and if a Ontology node does not describe a protein, then all its descendant Ontology nodes must not describe it”. The TPG constraint governs the structure of the Ontology-DAG and the inference process in the associated Forney Factor Graph.
}

\value{
A binary matrix with \emph{(n-1)^2+1} rows by \emph{n+1} columns where \emph{n = attr}
}
\references{
Tanoue J, Yoshikawa M, Uemura S (2012). The Gene Around GO viewer. Bioinformatics 18(12): 1705–1706.

Spetale FE, Tapia E, Krsticevic F, Roda F, Bulacio P (2016). A Factor Graph Approach to Automated GO Annotation. PLOS ONE 11(1): e0146986
}
\author{
Flavio E. Spetale, Pilar Bulacio and Javier Murillo <spetale@cifasis-conicet.gov.ar>
}

\examples{
tableTPG(3)
}

\keyword{TPG}
