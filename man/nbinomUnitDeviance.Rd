\name{nbinomUnitDeviance}
\alias{nbinomUnitDeviance}

\title{Negative Binomial Unit Deviance}

\description{
Compute unit deviances for the negative binomial distribution.
}

\usage{
nbinomUnitDeviance(y, mean, dispersion = 0)
}

\arguments{
\item{y}{numeric vector or matrix containing negative binomial counts.}

\item{mean}{numeric vector or matrix of expected values. If a matrix, then of same dimensions as \code{y}.}

\item{dispersion}{numeric vector or matrix of negative binomial dispersions.
Can be a scalar, a vector of length \code{nrow(y)} or a matrix of same dimensions as \code{y}.}
}

\details{
The unit deviance of the negative binomial distribution is a measure of the distance between \code{y} and {mean}.
If \code{mean} and \code{dispersion} are the true mean and dispersion of the negative binomial distribution, then the unit deviance follows an approximate chisquare distribution on 1 degree of freedom.

This function computes the unit deviance for each \code{y} observation.
Care is taken to ensure accurate computation in limiting cases when the dispersion is near zero or \code{mean*dispersion} is very large.
}

\value{
Numeric vector or matrix of the same size as \code{y} containing unit deviances.
}

\references{
Dunn PK, Smyth GK (2018).
\emph{Generalized linear models with examples in R}.
Springer, New York, NY.
\doi{10.1007/978-1-4419-0118-7}
ISBN: 978-1-4419-0118-7.

Jorgensen B (2013).
Generalized linear models.
\emph{Encyclopedia of Environmetrics} 3, Wiley.
\doi{10.1002/9780470057339.vag010.pub2}

McCarthy DJ, Chen Y, Smyth GK (2012).
Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
\emph{Nucleic Acids Research} 40, 4288-4297.
\doi{10.1093/nar/gks042}
}

\author{Gordon Smyth, Yunshun Chen, Aaron Lun.  C++ code by Aaron Lun.}

\examples{
y <- 1:4
names(y) <- letters[1:4]
nbinomUnitDeviance(y,mean=2.5,dispersion=0.2)
}

\concept{Model fit}