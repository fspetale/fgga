\name{CfData}
\alias{CfData}
\docType{data}
\title{
A set of characterized protein coding genes from the Cannis familiaris organism
annotated to a target GO subgraph considering both experimental and electronic
evidence.
}
\description{
The CfData dataset consists of a list containing the following:

\emph{$dxCf}: characterizations of \emph{6962} protein coding genes in terms of \emph{72} physico-chemical properties of their amino acid sequences. These sequences, obtained from the Uniprot database, are annotated to \emph{36} GO-terms of the GO Molecular Function (GO-MF) ontology subdomain.

\emph{$tableCfGO}: a set of \emph{6962} protein coding genes annotated to GO-MF target classes. Genes are identified by their \bold{Uniprot} ID mappings which are obtained with the org.Cf.eg.db annotation package set to work with both experimental and electronic evidence. Additionally, only those GO-MF terms with at least \emph{500} annotated genes were preserved.

\emph{$graphCfGO}: the target \bold{GO-MF} subgraph obtained with the org.Cf.eg.db annotation package set to work with the set of \bold{GO-MF} target classes.

\emph{$indexGO}: two arrays of Uniprot ID mappings defining the train-test partition of the set \emph{6962} protein coding genes annotated to \bold{GO-MF} terms.

\emph{$nodesGO}: labels of the GO-MF subgraph.

\emph{$varianceGOs}: a vector labeled with the variance of each \bold{GO-MF} term.
}
\usage{data("CfData")}
\format{
    A list with five named entries containg:
    \describe{
    \item{dxCf}{A matrix (6962 rows x 72 columns) containing the characterized proteins.}
    \item{graphCfGO}{An adjacency binary matrix (36 rows x 36 columns) corresponding to the GO-MF subgraph.}
    \item{indexGO}{A list with two named entries: indexTrain and indexTest each containing a numeric vector.}
    \item{tableCfGO}{A binary matrix (6962 rows x 36 columns) containing GOs associated with a protein.}
    \item{nodesGO}{A numerical vector containing the nodes of the GO-MF subgraph.}
    }
}


\source{
Uniprot Taxonomy: 9615

\url{https://www.uniprot.org/uniprot/?query=taxonomy:9615}

Package: org.Cf.eg.db - Version: 3.8.2

\url{https://bioconductor.org/packages/org.Cf.eg.db/}
}

\examples{
data(CfData)

## list objects included
ls(CfData)
# [1] "dxCf"  "graphCfGO" "indexGO"   "nodesGO"   "tableCfGO"

# Physico-chemical properties of each protein
head(CfData[["dxCf"]])

# GO-MF node labels, GO-terms, of each protein
head(CfData[["tableCfGO"]])

}
\keyword{datasets}
