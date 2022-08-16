\name{createFolds}
\alias{createFolds}
\title{
Data splitting function useful for binary classification tasks
}
\description{
createFolds splits binary classification data into k-folds.
}

\usage{
createFolds(target, k_fold = 10)
}

\arguments{
    \item{target}{A binary vector of a Ontology class}

    \item{k_fold}{An integer for the number of folds}
}
\details{
A random sampling is performed on binary classification data. A set of k data folds reflecting the original class balance is obtained.
}
\value{
list of row position integers corresponding to the training data
}
\references{
Hyndman and Athanasopoulos (2013), Forecasting: principles and practice. \url{https://www.otexts.org/fpp}
}
\author{
Flavio E. Spetale and Pilar Bulacio <spetale@cifasis-conicet.gov.ar>
}

\examples{
data(CfData)

createFolds(CfData[["tableCfGO"]][ ,"GO:0005515"], k_fold = 2)}

