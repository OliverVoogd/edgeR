\name{predFC}
\alias{predFC}
\alias{predFC.DGEList}
\alias{predFC.SummarizedExperiment}
\alias{predFC.default}

\title{Predictive log-fold changes}

\description{Computes estimated coefficients for a NB glm in such a way that the log-fold-changes are shrunk towards zero.}

\usage{
\method{predFC}{DGEList}(y, design, prior.count=0.125, offset=NULL, dispersion=NULL, 
       weights=NULL, \dots)
\method{predFC}{SummarizedExperiment}(y, design, prior.count=0.125, offset=NULL, dispersion=NULL, 
       weights=NULL, \dots)
\method{predFC}{default}(y, design, prior.count=0.125, offset=NULL, dispersion=0, 
       weights=NULL, \dots)
}

\arguments{
\item{y}{a matrix of counts, or a \code{DGEList} object, or a \code{SummarizedExperiment} object} 

\item{design}{the design matrix for the experiment}

\item{prior.count}{the average prior count to be added to each observation. Larger values produce more shrinkage.}

\item{offset}{numeric vector or matrix giving the offset in the log-linear model predictor, as for \code{\link{glmFit}}. Usually equal to log library sizes.}

\item{dispersion}{numeric vector of negative binomial dispersions.}

\item{weights}{optional numeric matrix giving observation weights}

\item{\ldots}{other arguments are passed to \code{glmFit}.}
}

\details{
This function computes predictive log-fold changes (pfc) for a NB GLM.
The pfc are posterior Bayesian estimators of the true log-fold-changes.
They are predictive of values that might be replicated in a future experiment.

Specifically, the function adds a small prior count to each observation before fitting the GLM (see \code{\link{addPriorCount}} for details).
The actual prior count that is added is proportion to the library size.
This has the effect that any log-fold-change that was zero prior to augmentation remains zero and non-zero log-fold-changes are shrunk towards zero.

The prior counts can be viewed as equivalent to a prior belief that the log-fold changes are small, and the output can be viewed as posterior log-fold-changes from this Bayesian viewpoint.
The output coefficients are called \emph{predictive} log fold-changes because, depending on the prior, they may be a better prediction of the true log fold-changes than the raw estimates.

Log-fold changes for genes with low counts are shrunk more than those for genes with high counts.
In particular, infinite log-fold-changes arising from zero counts are avoided.
The exact degree to which this is done depends on the negative binomail dispersion.
}

\value{
Numeric matrix of (shrunk) linear model coefficients on the log2 scale.
}

\author{Belinda Phipson and Gordon Smyth}

\references{
Phipson, B. (2013).
\emph{Empirical Bayes modelling of expression profiles and their associations}.
PhD Thesis. University of Melbourne, Australia.
\url{http://repository.unimelb.edu.au/10187/17614}
}

\seealso{
\code{\link{glmFit}}, \code{\link{exactTest}},
\code{\link{addPriorCount}}
}

\examples{
# generate counts for a two group experiment with n=2 in each group and 100 genes
disp <- 0.1
y <- matrix(rnbinom(400,size=1/disp,mu=4), nrow=100, ncol=4)
y <- DGEList(y, group=c(1,1,2,2))
design <- model.matrix(~group, data=y$samples)

#estimate the predictive log fold changes
predlfc <- predFC(y, design, dispersion=disp, prior.count=1)
logfc <- predFC(y,design,dispersion=disp, prior.count=0)
logfc.truncated <- pmax(pmin(logfc,100),-100)

#plot predFC's vs logFC's
plot(predlfc[,2], logfc.truncated[,2],
     xlab="Predictive log fold changes", ylab="Raw log fold changes")
abline(a=0,b=1)
}

\concept{Model fit}
