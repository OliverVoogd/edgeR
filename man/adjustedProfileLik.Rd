\name{adjustedProfileLik}
\alias{adjustedProfileLik}

\title{Adjusted Profile Likelihood for the Negative Binomial Dispersion Parameter}

\description{
Compute adjusted profile log-likelihoods for the dispersion parameters of genewise negative binomial glms.
}

\usage{
adjustedProfileLik(dispersion, y, design, offset, weights=NULL, adjust=TRUE, 
            start=NULL, get.coef=FALSE)
}

\arguments{
\item{dispersion}{numeric scalar or vector of dispersions.}
\item{y}{numeric matrix of counts.}
\item{design}{numeric matrix giving the design matrix.}
\item{offset}{numeric matrix of same size as \code{y} giving offsets for the log-linear models.  Can be a scalor or a vector of length \code{ncol(y)}, in which case it is expanded out to a matrix.}
\item{weights}{optional numeric matrix giving observation weights.}
\item{adjust}{logical, if \code{TRUE} then Cox-Reid adjustment is made to the log-likelihood, if \code{FALSE} then the log-likelihood is returned without adjustment.}
\item{start}{numeric matrix of starting values for the GLM coefficients, to be passed to \code{\link{glmFit}}.}
\item{get.coef}{logical, specifying whether fitted GLM coefficients should be returned.}
}

\value{
If \code{get.coef==FALSE}, a vector of adjusted profile log-likelihood values is returned containing one element for each row of \code{y}.

Otherwise, a list is returned containing \code{apl}, the aforementioned vector of adjusted profile likelihoods, and \code{beta}, the numeric matrix of fitted GLM coefficients.
}

\details{
For each row of data, compute the adjusted profile log-likelihood for the dispersion parameter of the negative binomial glm.
The adjusted profile likelihood is described by McCarthy et al (2012) and is based on the method of Cox and Reid (1987).

The adjusted profile likelihood is an approximation to the log-likelihood function, conditional on the estimated values of the coefficients in the NB log-linear models.
The conditional likelihood approach is a technique for adjusting the likelihood function to allow for the fact that nuisance parameters have to be estimated in order to evaluate the likelihood.
When estimating the dispersion, the nuisance parameters are the coefficients in the log-linear model.

This implementation calls the LAPACK library to perform the Cholesky decomposition during adjustment estimation.

The purpose of \code{start} and \code{get.coef} is to allow hot-starting for multiple calls to \code{adjustedProfileLik}, when only the \code{dispersion} is altered.
Specifically, the returned GLM coefficients from one call with \code{get.coef==TRUE} can be used as the \code{start} values for the next call.

The \code{weights} argument is interpreted in terms of averages.
Each value of \code{y} is assumed to be the average of \code{n} independent and identically distributed NB counts, where \code{n} is given by the weight.
This assumption can generalized to fractional weights.
}

\references{
Cox, DR, and Reid, N (1987). Parameter orthogonality and approximate conditional inference. \emph{Journal of the Royal Statistical Society Series B} 49, 1-39.

McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
\emph{Nucleic Acids Research} 40, 4288-4297.
\doi{10.1093/nar/gks042}
}

\author{Yunshun Chen, Gordon Smyth, Aaron Lun}
\examples{
y <- matrix(rnbinom(30, mu=10, size=20), 10, 3)
design <- matrix(1, 3, 1)
dispersion <- 0.05
adjustedProfileLik(dispersion, y, design, offset=0)
}

\seealso{
\code{\link{glmFit}}
}

\concept{Dispersion estimation}
