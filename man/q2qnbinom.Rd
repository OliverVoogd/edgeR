\name{q2qnbinom}
\alias{q2qpois}
\alias{q2qnbinom}
\title{Quantile to Quantile Mapping between Negative-Binomial Distributions}

\description{Interpolated quantile to quantile mapping between negative-binomial distributions with the same dispersion but different means.
The Poisson distribution is a special case.}

\usage{
q2qpois(x, input.mean, output.mean)
q2qnbinom(x, input.mean, output.mean, dispersion=0)
}

\arguments{
\item{x}{numeric matrix of counts.}
\item{input.mean}{numeric matrix of population means for \code{x}.  If a vector, then of the same length as \code{nrow(x)}.}
\item{output.mean}{numeric matrix of population means for the output values.  If a vector, then of the same length as \code{nrow(x)}.}
\item{dispersion}{numeric scalar, vector or matrix giving negative binomial dispersion values.}
}

\details{
This function finds the quantile with the same left and right tail probabilities relative to the output mean as \code{x} has relative to the input mean.
\code{q2qpois} is equivalent to \code{q2qnbinom} with \code{dispersion=0}.

In principle, \code{q2qnbinom} gives similar results to calling \code{pnbinom} followed by \code{qnbinom} as in the example below.
However this function avoids infinite values arising from rounding errors and does appropriate interpolation to return continuous values.

\code{q2qnbinom} is called by \code{\link{equalizeLibSizes}} to perform quantile-to-quantile normalization.
}

\value{numeric matrix of same dimensions as \code{x}, with \code{output.mean} as the new nominal population mean.}

\seealso{
\code{\link{equalizeLibSizes}}
}

\author{Gordon Smyth}

\examples{
x <- 15
input.mean <- 10
output.mean <- 20
dispersion <- 0.1
q2qnbinom(x,input.mean,output.mean,dispersion)

# Similar in principle:
qnbinom(pnbinom(x,mu=input.mean,size=1/dispersion),mu=output.mean,size=1/dispersion)
}

\concept{Normalization}
