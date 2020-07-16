\name{estimateGLMRobustDisp}
\alias{estimateGLMRobustDisp}

\title{Empirical Robust Bayes Tagwise Dispersions for Negative Binomial GLMs using Observation Weights}

\description{
Compute a robust estimate of the negative binomial dispersion parameter for each gene, with expression levels specified by a log-linear model, using observation weights.  These observation weights will be stored and used later for estimating regression parameters.
}

\usage{
estimateGLMRobustDisp(y, design = NULL, prior.df = 10, update.trend = TRUE,
                      trend.method = "bin.loess", maxit = 6, k = 1.345,
                      residual.type = "pearson", verbose = FALSE,
                      record = FALSE)
}

\arguments{
\item{y}{a \code{DGEList} object.}
\item{design}{numeric design matrix, as for \code{\link{glmFit}}.}
\item{prior.df}{prior degrees of freedom.}
\item{update.trend}{logical. Should the trended dispersion be re-estimated at each iteration?} 
\item{trend.method}{method (low-level function) used to estimated the trended dispersions. \code{\link{estimateGLMTrendedDisp}}}
\item{maxit}{maximum number of iterations for weighted \code{\link{estimateGLMTagwiseDisp}}.}
\item{k}{the tuning constant for Huber estimator. If the absolute value of residual (r) is less than k, its observation weight is 1, otherwise \code{k/abs(r)}.} 
\item{residual.type}{type of residual (r) used for estimation observation weight}
\item{verbose}{logical.  Should verbose comments be printed?}
\item{record}{logical.  Should information for each iteration be recorded (and returned as a list)?}
}

\value{
\code{estimateGLMRobustDisp} produces a \code{DGEList} object, which contains the (robust) genewise dispersion parameter estimate for each gene for the negative binomial model that maximizes the weighted Cox-Reid adjusted profile likelihood, as well as the observation weights.  The observation weights are calculated using residuals and the Huber function.

Note that when \code{record=TRUE}, a simple list of \code{DGEList} objects is returned, one for each iteration (this is for debugging or tracking purposes).
}

\details{
Moderation of dispersion estimates towards a trend can be sensitive to outliers, resulting in an increase in false positives. That is, since the dispersion estimates are moderated downwards toward the trend and because the regression parameter estimates may be affected by the outliers, some genes are incorrectly deemed to be significantly differentially expressed. This function uses an iterative procedure where weights are calculated from residuals and estimates are made after re-weighting.

The robustly computed genewise estimates are reported in the \code{tagwise.dispersion} vector of the returned \code{DGEList}.
The terms `tag' and `gene' are synonymous in this context.

Note: it is not necessary to first calculate the common, trended and genewise dispersion estimates.  If these are not available, the function will first calculate this (in an unweighted) fashion.
}

\references{
Zhou X, Lindsay H, Robinson MD (2014). Robustly detecting differential expression in RNA sequencing data using observation weights. Nucleic Acids Research, 42(11), e91.
}

\author{Xiaobei Zhou, Mark D. Robinson}
\examples{
y <- matrix(rnbinom(100*6,mu=10,size=1/0.1),ncol=6)
d <- DGEList(counts=y,group=c(1,1,1,2,2,2),lib.size=c(1000:1005))
d <- calcNormFactors(d)
design <- model.matrix(~group, data=d$samples) # Define the design matrix for the full model
d <- estimateGLMRobustDisp(d, design)
summary(d$tagwise.dispersion)
}

\seealso{
This function calls
\code{\link{estimateGLMTrendedDisp}}
and
\code{\link{estimateGLMTagwiseDisp}}.
}

\concept{Dispersion estimation}
