\name{sumTechReps}
\alias{sumTechReps}
\alias{sumTechReps.default}
\alias{sumTechReps.DGEList}
\alias{sumTechReps.SummarizedExperiment}

\title{Sum Over Replicate Samples}
\description{
Condense the columns of a matrix or DGEList object so that counts are summed over technical replicate samples.
}
\usage{
\method{sumTechReps}{default}(x, ID=colnames(x), \dots)
\method{sumTechReps}{DGEList}(x, ID=colnames(x), \dots)
\method{sumTechReps}{SummarizedExperiment}(x, ID, \dots)
}
\arguments{
  \item{x}{a numeric matrix or \code{DGEList} object.}
  \item{ID}{sample identifier.}
  \item{\dots}{other arguments are not currently used.}
}
\details{
A new matrix or \code{DGEList} object is computed in which the counts for technical replicate samples are replaced by their sums.
}
\value{
A data object of the same class as \code{x} with a column for each unique value of \code{ID}.
For a \code{SummarizedExperiment} object, it is converted into a \code{DGEList} object.
Columns are in the same order as the ID values first occur in the ID vector.
}
\author{Gordon Smyth and Yifang Hu}
\seealso{
   \code{\link{rowsum}}.
}

\examples{
x <- matrix(rpois(8*3,lambda=5),8,3)
colnames(x) <- c("a","a","b")
sumTechReps(x)
}
