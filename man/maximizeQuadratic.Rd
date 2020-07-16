\name{maximizeQuadratic}
\alias{maximizeQuadratic}

\title{Maximize a function given a table of values by quadratic interpolation.}

\description{
Maximize a function given a table of values by quadratic interpolation.
}

\usage{
maximizeQuadratic(y, x=1:ncol(y))
}

\arguments{ 
\item{y}{numeric matrix of response values.}

\item{x}{numeric matrix of inputs of the function of same dimension as \code{y}.  If a vector, must be a row vector of length equal to \code{ncol(y)}.}
}

\details{
For each row of \code{y}, finds the three \code{x} values bracketing the maximum of \code{y}, interpolates a quadatric
polyonomial through these \code{y} for these three values and solves for the location of the maximum of the polynomial.
}

\value{
numeric vector of length equal to \code{nrow(y)} giving the x-value at which \code{y} is maximized.
}

\author{Yunshun Chen and Gordon Smyth}

\examples{
y <- matrix(rnorm(5*9),5,9)
maximizeQuadratic(y)
}

\seealso{
\code{\link{maximizeInterpolant}}
}
