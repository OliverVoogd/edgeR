\name{cutWithMinN}
\alias{cutWithMinN}

\title{Cut numeric vector into non-empty intervals}

\description{
Discretizes a numeric vector.
Divides the range of \code{x} into intervals, so that each interval contains a minimum number of values, and codes the values in \code{x} according to which interval they fall into.
}

\usage{cutWithMinN(x, intervals=2, min.n=1)}

\arguments{ 
	\item{x}{numeric vector.}
	\item{intervals}{number of intervals required.}
	\item{min.n}{minimum number of values in any interval. Must be less than or equal to \code{length(x)/intervals}.}
}

\value{
A list with components:
\item{group}{integer vector of same length as \code{x} indicating which interval each value belongs to.}
\item{breaks}{numeric vector of length \code{intervals+1} giving the left and right limits of each interval.}
}

\details{
This function strikes a compromise between the base functions \code{cut},
which by default cuts a vector into equal length intervals,
and \code{quantile}, which is suited to finding equally populated intervals.
It finds a partition of the \code{x} values that is as close as possible to equal length intervals while keeping at least \code{min.n} values in each interval.

Tied values of \code{x} are broken by random jittering, so the partition may vary slightly from run to run if there are many tied values.
}

\author{Gordon Smyth}

\seealso{
\code{\link{cut}}, \code{\link{quantile}}.
}

\examples{
x <- c(1,2,3,4,5,6,7,100)
cutWithMinN(x,intervals=3,min.n=1)
}

