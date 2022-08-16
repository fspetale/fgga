\name{fMeasure}
\alias{fMeasures}
\alias{fHierarchicalMeasures}
\alias{fMeasuresByLevel}

\title{
Individual and hierarchical F-measures
}
\description{
Set of functions to compute the individual and hierarchical F-score, precision, recall.
}
\usage{
    fMeasures(target, predicted, cutoff = 0.5)
    fMeasuresByLevel(target, predicted, graphOnto, cutoff = 0.5)
    fHierarchicalMeasures(target, predicted, graphOnto, cutoff = 0.5)
}

\arguments{
    \item{target}{
    A binary matrix with `n' proteins (rows) by `m' Ontology node labels (columns) corresponding to the target of ontology terms where 0 stands for negative and 1 for positive.
    }
    \item{predicted}{
    A real matrix with `n' proteins (rows) by `m' Ontology node labels (columns) corresponding to the predicted terms.
    }
    \item{graphOnto}{
    A graphNEL graph with `m' Ontology node labels.
    }
    \item{cutoff}{
    A real value to divide the predicted terms into positive and negative. The predicted values higher than the cutoff will be taken as positive.
    }
}

\details{
    \code{fMeasures} computes the F-score, precision, recall, specificity and accuracy for each ontological term.

    \code{fMeasuresByLevel} computes F-score, precision, recall, specificity and accuracy for all ontological terms belongs to graph. The levels are calculated as the maximum distance between two terms of the graph.

    \code{fHierarchicalMeasures} computes the hierarchical F-score, precision, recall for the predicted terms of a set of proteins.
}

\value{
    \code{fMeasures} and \code{fMeasuresByLevel} returns a list of two elements where the first element is a named vector with six attributes while the second element is an array of 'm' ontological terms by six attributes. The 6 attributes are:
    \item{Prec:}{Precision}
    \item{Recall:}{Recall}
    \item{Specif:}{Specificity}
    \item{Fmeasure:}{F-score}
    \item{Acc:}{Accuracy}
    \item{nPositive:}{Number of positive samples}

    \code{fHierarchicalMeasures} returns a list of five elements:
    \item{HP:}{Hierarchical Precision}
    \item{HR:}{Hierarchical Recall}
    \item{HF:}{Hierarchical F-score}
    \item{nSample: }{Number of proteins evaluated}
    \item{noEvalSample: }{Named vector of proteins not evaluated}
}

\references{
Verspoor K, Cohn J, Mnizewski S, C J. A categorization approach to automated ontological function annotation. Protein Science. 2006;15:1544–1549.
}

\author{
Flavio E. Spetale <spetale@cifasis-conicet.gov.ar>
}

\examples{
data(CfData)

predGO <- matrix(runif(360, 0, 1),10,36, dimnames=list(rownames(
        CfData[["tableCfGO"]])[seq_len(10)], colnames(CfData[["tableCfGO"]])))

fMeasures(CfData[["tableCfGO"]][seq_len(10), ], predGO, cutoff = 0.5)

mygraphGO <- as(CfData[["graphCfGO"]], "graphNEL")

fHierarchicalMeasures(CfData[["tableCfGO"]][seq_len(10), ], predGO, mygraphGO,
cutoff = 0.5)
}
