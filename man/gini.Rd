\name{gini}
\alias{gini}
\title{Gini dispersion index}
\description{
Gini index for each column of a matrix.
}
\usage{
gini(x)
}
\arguments{
  \item{x}{a non-negative numeric matrix, or an object that can be coerced to such a matrix by \code{as.matrix}.}
}
\details{
The Gini coefficient or index is a measure of inequality or diversity.
It is zero if all the values of \code{x} are equal.
It reaches a maximum value of \code{1/nrow(x)} when all values are zero except for one.

The Gini index is only interpretable for non-negative quantities.
It is not meaningful if \code{x} contains negative values.
}
\value{
Numeric vector of length \code{ncol(x)}.
}
\references{
\url{https://en.wikipedia.org/wiki/Gini_coefficient}.
}
\examples{
x <- matrix(rpois(20,lambda=5),10,2)
gini(x)
}
\author{Gordon Smyth}

\concept{Data exploration}
