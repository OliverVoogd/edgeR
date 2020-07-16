\name{thinCounts}
\alias{thinCounts}

\title{Binomial or Multinomial Thinning of Counts}

\description{
Reduce the size of Poisson-like counts by binomial thinning.
}

\usage{thinCounts(x, prob=NULL, target.size=min(colSums(x)))}

\arguments{ 
	\item{x}{numeric vector or array of non-negative integers.}
	\item{prob}{numeric scalar or vector of same length as \code{x}, the expected proportion of the events to keep.}
	\item{target.size}{integer scale or vector of the same length as \code{NCOL{x}}, the desired total column counts.
	Must be not greater than column sum of \code{x}.	Ignored if \code{prob} is not \code{NULL}.}
}

\value{
A vector or array of the same dimensions as \code{x}, with thinned counts.
}

\details{
If \code{prob} is not \code{NULL}, then this function calls \code{rbinom} with \code{size=x} and \code{prob=prob} to generate the new counts.
This is classic binomial thinning.
The new column sums are random, with expected values determined by \code{prob}.

If \code{prob} is \code{NULL}, then this function does multinomial thinning of the counts to achieve specified column totals.
The default behavior is to thin the columns to have the same column sum, equal to the smallest column sum of \code{x}.

If the elements of \code{x} are Poisson, then binomial thinning produces new Poisson random variables with expected values reduced by factor \code{prob}.
If the elements of each column of \code{x} are multinomial, then multinomial thinning produces a new multinomial observation with a reduced sum.
}

\author{Gordon Smyth}

\examples{
x <- rpois(10,lambda=10)
thinCounts(x,prob=0.5)
}

\concept{Data exploration}


