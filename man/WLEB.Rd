\name{WLEB}
\alias{WLEB}

\title{Weighted Likelihood Empirical Bayes}

\description{
Compute empirical Bayes moderated parameter estimators using a weighted likelihood approach.
}

\usage{
WLEB(theta, loglik, prior.n = 5, covariate = NULL, trend.method = "locfit", span = NULL, 
     overall = TRUE, trend = TRUE, individual = TRUE, m0 = NULL, m0.out = FALSE)
}

\arguments{
\item{theta}{numeric vector of values of the parameter at which the log-likelihoods are calculated.}
\item{loglik}{numeric matrix of log-likelihood of all the candidates at those values of parameter.}
\item{prior.n}{numeric scaler, estimate of the prior weight, i.e. the smoothing parameter that indicates the weight to put on the common likelihood compared to the individual's likelihood.}
\item{covariate}{numeric vector of values across which a parameter trend is fitted}
\item{trend.method}{method for estimating the parameter trend. Possible values are \code{"none"}, \code{"movingave"}, \code{"loess"}, \code{"locfit"} or \code{"locfit.mixed"}. The latter options cause \code{movingAverageByCol()}, \code{loessBycol()} and \code{locfitByCol()} respectively to be called to smooth \code{loglik} by the covariate values. The \code{"locfit.mixed"} method is the same as \code{"locfit"} but uses a polynomial of degree 1 for lowly expressed genes.}
\item{span}{width of the smoothing window, as a proportion of the data set.}
\item{overall}{logical, should a single value of the parameter which maximizes the sum of all the log-likelihoods be estimated?}
\item{trend}{logical, should a parameter trend (against the covariate) which maximizes the local shared log-likelihoods be estimated?}
\item{individual}{logical, should individual estimates of all the candidates after applying empirical Bayes method along the trend be estimated?}
\item{m0}{numeric matrix of local shared log-likelihoods. If \code{NULL}, it will be calculated using the method selected by \code{trend.method}.}
\item{m0.out}{logical, should local shared log-likelihoods be included in the output?}
}

\details{
This function implements a very general empirical Bayes strategy outlined by McCarthy et al (2012).
McCarthy et el used the method to estimate genewise negative binomial dispersion parameters, but here the method is generalized to apply to any parameter and any likelihood function.
The method gives similar results to parametric empirical Bayes with a conjugate prior, when a conjugate prior exists, but does not require the prior distribution to be specified.
The prior distribution is instead inferred from the pooled likelihood, i.e., from the likelihood that would be arise from pooling all the cases or from pooling all the cases with similar covariate values.

The function assumes a series of cases.
Each case leads to data set from which a parameter (theta) is to be estimated.
For each case, the log-likelihood function has been evaluated over a grid of possible values for theta.
The function takes as input the matrix of log-likelihood values where the rows correspond to cases and the columns correspond to putative parameter values.

Each case is associated with a covariate value that might affect theta.
The "overall" parameter estimate is the maximum likelihood estimator of theta that arises if the likelihood is averaged over cases and then maximized over theta.
The "trend" parameter estimates are estimates for theta that arise if each column of \code{loglik} is replaced by a smooth trend with respect to the covariate.
The "individual" parameter estimate for each case is a compromise between the maximum likelihood estimate for that case alone and a global parameter estimate computed from all the cases, the latter being either the overall estimate (if \code{trend.method="none"} or the trend estimate (otherwise).
}

\value{
A list with the following:
	\item{overall}{the parameter estimate that maximizes the sum of all the log-likelihoods.}
	\item{trend}{the estimated trended parameters against the covariate.}
	\item{individual}{the individual estimates of all the candidates after applying empirical Bayes method along the trend.}
	\item{shared.loglik}{the estimated numeric matrix of local shared log-likelihoods}
}

\author{Yunshun Chen, Gordon Smyth}

\references{
McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
\emph{Nucleic Acids Research} 40, 4288-4297.
\doi{10.1093/nar/gks042}
}

\examples{
y <- matrix(rpois(100, lambda=10), ncol=4)
theta <- 7:14
loglik <- matrix(0,nrow=nrow(y),ncol=length(theta))
for(i in 1:nrow(y))
for(j in 1:length(theta))
	loglik[i,j] <- sum(dpois(y[i,], theta[j] ,log=TRUE))
covariate <- log(rowSums(y))
out <- WLEB(theta, loglik, prior.n=3, covariate)
out
}

\seealso{
\code{\link{locfitByCol}}, \code{\link{movingAverageByCol}} and \code{\link{loessByCol}} implement the local fit, moving average or loess smoothers.
}

\concept{Dispersion estimation}
