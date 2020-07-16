\name{nearestReftoX}
\alias{nearestReftoX}
\title{Find Nearest Element of Reference for each Element of X}

\description{
Find nearest element of a sorted reference vector and to each element of x.
}

\usage{
nearestReftoX(x, reference, \ldots)
}

\arguments{
  \item{x}{numeric vector.}
  \item{reference}{numeric vector, sorted in increasing order.}
  \item{\dots}{other arguments as passed to \code{findInterval}.}
}

\details{
This function finds the element of a reference table (\code{reference}) that is closest to each element of an incoming vector (\code{x}).

The function is a simple wrapper for \code{findInterval} in the base package.
It calls \code{findInterval} with \code{vec} equal to the mid-points between the reference values.
}

\value{
Integer vector giving indices of elements of \code{reference}.
}

\author{Gordon Smyth}

\seealso{
\code{\link{findInterval}}
}

\examples{
nearestReftoX(c(-10,0.5,0.6,2,3), reference = c(-1,0,2))
}
