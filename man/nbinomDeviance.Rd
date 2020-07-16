\name{nbinomDeviance}
\alias{nbinomDeviance}

\title{Negative Binomial Deviance}

\description{
Fit the same log-link negative binomial or Poisson generalized linear model (GLM) to each row of a matrix of counts.
}

\usage{
nbinomDeviance(y, mean, dispersion=0, weights=NULL)
}

\arguments{
\item{y}{numeric matrix containing the negative binomial counts, with rows for genes and columns for libraries. A vector will be treated as a matrix with one row.}

\item{mean}{numeric matrix of expected values, of same dimension as \code{y}. A vector will be treated as a matrix with one row.}

\item{dispersion}{numeric vector or matrix of negative binomial dispersions, as for \code{\link{glmFit}}.
Can be a scalar, a vector of length equal to \code{nrow(y)}, or a matrix of same dimensions as \code{y}.}

\item{weights}{numeric vector or matrix of non-negative weights, as for \code{\link{glmFit}}.
Can be a scalar, a vector of length equal to \code{ncol(y)}, or a matrix of same dimensions as \code{y}.}
}

\details{
Computes the total residual deviance for each row of \code{y}, i.e., weighted row sums of the unit deviances.

Care is taken to ensure accurate computation in limiting cases when the dispersion is near zero or \code{mean*dispersion} is very large.
}

\value{
\code{nbinomDeviance} returns a numeric vector of length equal to the number of rows of \code{y}.
}

\references{
Jorgensen, B. (2013).
Generalized linear models. Encyclopedia of Environmetrics 3, Wiley.
\url{http://onlinelibrary.wiley.com/doi/10.1002/9780470057339.vag010.pub2/abstract}.

McCarthy, DJ, Chen, Y, Smyth, GK (2012).
Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
\emph{Nucleic Acids Research} 40, 4288-4297.
\url{https://doi.org/10.1093/nar/gks042}
}

\author{Gordon Smyth, Yunshun Chen, Aaron Lun.  C++ code by Aaron Lun.}

\seealso{
\code{\link{nbinomUnitDeviance}}
}

\examples{
y <- matrix(1:6,3,2)
mu <- matrix(3,3,2)
nbinomDeviance(y,mu,dispersion=0.2)
}

\concept{Model fit}
