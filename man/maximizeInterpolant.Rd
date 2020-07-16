\name{maximizeInterpolant}
\alias{maximizeInterpolant}

\title{Maximize a function given a table of values by spline interpolation.}

\description{
Maximize a function given a table of values by spline interpolation.
}

\usage{
maximizeInterpolant(x, y)
}

\arguments{
\item{x}{numeric vector of the inputs of the function.}

\item{y}{numeric matrix of function values at the values of \code{x}.
Columns correspond to \code{x} values and each row corresponds to a different function to be maximized.}
}

\value{
numeric vector of input values at which the function maximums occur.
}

\details{
Calculates the cubic spline interpolant for each row the method of Forsythe et al (1977) using the function \code{fmm_spline} from \code{splines.c} in the \code{stats} package).
Then calculates the derivatives of the spline segments adjacant to the input with the maximum function value.
This allows identification of the maximum of the interpolating spline.
}

\author{Aaron Lun, improving on earlier code by Gordon Smyth}
\examples{
x <- seq(0,1,length=10)
y <- rnorm(10,1,1)
maximizeInterpolant(x,y)
}

\references{
Forsythe, G. E., Malcolm, M. A. and Moler, C. B. (1977).
\emph{Computer Methods for Mathematical Computations}, Prentice-Hall.
}
