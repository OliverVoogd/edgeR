\name{dropEmptyLevels}
\alias{dropEmptyLevels}
\title{Drop Levels of a Factor that Never Occur}
\description{
Reform a factor so that only necessary levels are kept.
}
\usage{
dropEmptyLevels(x)
}
\arguments{
  \item{x}{a factor or a vector to be converted to a factor.}
}
\details{
In general, the levels of a factor, \code{levels(x)}, may include values that never actually occur.
This function drops any levels of that do not occur.

If \code{x} is not a factor, then the function returns \code{factor(x)}.
If \code{x} is a factor, then the function returns the same value as \code{factor(x)} or \code{x[,drop=TRUE]} but somewhat more efficiently.
}
\value{
A factor with the same values as \code{x} but with a possibly reduced set of levels.
}
\author{Gordon Smyth}
\seealso{
   \code{\link{factor}}.
}

\examples{
x <- factor(c("a","b"), levels=c("c","b","a"))
x
dropEmptyLevels(x)
}
