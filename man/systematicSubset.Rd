\name{systematicSubset}
\alias{systematicSubset}

\title{Take a systematic subset of indices.}

\description{
Take a systematic subset of indices stratified by a ranking variable.
}

\usage{
systematicSubset(n, order.by)
}

\arguments{ 
\item{n}{integer giving the size of the subset.}

\item{order.by}{numeric vector of the values by which the indices are ordered.}
}

\value{
\code{systematicSubset} returns a vector of size \code{n}.
}


\author{Gordon Smyth}

\examples{
y <- rnorm(100, 1, 1)
systematicSubset(20, y)
}

\seealso{
\code{\link{order}}
}
