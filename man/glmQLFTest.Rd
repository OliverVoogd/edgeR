\name{glmQLFit}
\alias{glmQLFit}
\alias{glmQLFit.DGEList}
\alias{glmQLFit.SummarizedExperiment}
\alias{glmQLFit.default}
\alias{glmQLFTest}

\title{Genewise Negative Binomial Generalized Linear Models with Quasi-likelihood Tests}

\description{Fit a quasi-likelihood negative binomial generalized log-linear model to count data.
Conduct genewise statistical tests for a given coefficient or contrast.}

\usage{
\method{glmQLFit}{DGEList}(y, design=NULL, dispersion=NULL, abundance.trend=TRUE,
        robust=FALSE, winsor.tail.p=c(0.05, 0.1), \dots)
\method{glmQLFit}{SummarizedExperiment}(y, design=NULL, dispersion=NULL, abundance.trend=TRUE,
        robust=FALSE, winsor.tail.p=c(0.05, 0.1), \dots)
\method{glmQLFit}{default}(y, design=NULL, dispersion=NULL, offset=NULL, lib.size=NULL,
        weights=NULL, abundance.trend=TRUE, AveLogCPM=NULL, robust=FALSE,
        winsor.tail.p=c(0.05, 0.1), \dots)
glmQLFTest(glmfit, coef=ncol(glmfit$design), contrast=NULL, poisson.bound=TRUE)
}

\arguments{
\item{y}{a matrix of counts, or a \code{DGEList} object with (at least) elements \code{counts} (table of unadjusted counts) and \code{samples} (data frame containing information about experimental group, library size and normalization factor for the library size), or a \code{SummarizedExperiment} object with (at least) an element \code{counts} in its \code{assays}.}

\item{design}{numeric matrix giving the design matrix for the genewise linear models.}

\item{dispersion}{numeric scalar, vector or matrix of negative binomial dispersions. If \code{NULL}, then will be extracted from the \code{DGEList} object \code{y}, with order of precedence: trended dispersions, common dispersion, a constant value of 0.05.}

\item{offset}{numeric matrix of same size as \code{y} giving offsets for the log-linear models.  Can be a scalor or a vector of length \code{ncol(y)}, in which case it is expanded out to a matrix. If \code{NULL} will be computed by \code{getOffset(y)}.}

\item{lib.size}{numeric vector of length \code{ncol(y)} giving library sizes. Only used if \code{offset=NULL}, in which case \code{offset} is set to \code{log(lib.size)}. Defaults to \code{colSums(y)}.}

\item{weights}{numeric matrix of same size as \code{y} giving weights for the log-linear models. If \code{NULL}, will be set to unity for all observations.}

\item{abundance.trend}{logical, whether to allow an abundance-dependent trend when estimating the prior values for the quasi-likelihood multiplicative dispersion parameter.}

\item{AveLogCPM}{average log2-counts per million, the average taken over all libraries in \code{y}. If \code{NULL} will be computed by \code{aveLogCPM(y)}.}

\item{robust}{logical, whether to estimate the prior QL dispersion distribution robustly.}

\item{winsor.tail.p}{numeric vector of length 2 giving proportion to trim (Winsorize) from lower and upper tail of the distribution of genewise deviances when estimating the hyperparameters. Positive values produce robust empirical Bayes ignoring outlier small or large deviances. Only used when \code{robust=TRUE}.} 

\item{\dots}{other arguments are passed to \code{\link{glmFit}}.}

\item{glmfit}{a \code{DGEGLM} object, usually output from \code{glmQLFit}.}

\item{coef}{integer or character index vector indicating which coefficients of the linear model are to be tested equal to zero. Ignored if \code{contrast} is not \code{NULL}.}

\item{contrast}{numeric vector or matrix specifying one or more contrasts of the linear model coefficients to be tested equal to zero.}

\item{poisson.bound}{logical, if \code{TRUE} then the p-value returned will never be less than would be obtained for a likelihood ratio test with NB dispersion equal to zero.}
}

\details{
\code{glmQLFit} and \code{glmQLFTest} implement the quasi-likelihood (QL) methods of Lund et al (2012), with some enhancements and with slightly different glm, trend and FDR methods.
See Lun et al (2016) or Chen et al (2016) for tutorials describing the use of \code{glmQLFit} and \code{glmQLFit} as part of a complete analysis pipeline.
Another case study using \code{glmQLFit} and \code{glmQLFTest} is given in Section 4.7 of the edgeR User's Guide.

\code{glmQLFit} is similar to \code{glmFit} except that it also estimates QL dispersion values.
It calls the limma function \code{\link{squeezeVar}} to conduct empirical Bayes moderation of the genewise QL dispersions.
If \code{robust=TRUE}, then the robust hyperparameter estimation features of \code{squeezeVar} are used (Phipson et al, 2013).
If \code{abundance.trend=TRUE}, then a prior trend is estimated based on the average logCPMs.

\code{glmQLFit} gives special attention to handling of zero counts, and in particular to situations when fitted values of zero provide no useful residual degrees of freedom for estimating the QL dispersion (Lun and Smyth, 2017).
The usual residual degrees of freedom are returned as \code{df.residual} while the adjusted residual degrees of freedom are returned as \code{df.residuals.zeros}.

\code{glmQLFTest} is similar to \code{glmLRT} except that it replaces likelihood ratio tests with empirical Bayes quasi-likelihood F-tests.
The p-values from \code{glmQLFTest} are always greater than or equal to those that would be obtained from \code{glmLRT} using the same negative binomial dispersions.
}

\note{
The negative binomial dispersions \code{dispersion} supplied to \code{glmQLFit} and \code{glmQLFTest} must be based on a global model, that is, they must be either trended or common dispersions.
It is not correct to supply genewise dispersions because \code{glmQLFTest} estimates genewise variability using the QL dispersion.
}

\value{
\code{glmQLFit} produces an object of class \code{DGEGLM} with the same components as produced by \code{\link{glmFit}}, plus:
	\item{df.residual.zeros}{a numeric vector containing the number of effective residual degrees of freedom for each gene, taking into account any treatment groups with all zero counts.}
	\item{df.prior}{a numeric vector or scalar, giving the prior degrees of freedom for the QL dispersions.}
	\item{var.prior}{a numeric vector of scalar, giving the location of the prior distribution for the QL dispersions.}
	\item{var.post}{a numeric vector containing the posterior empirical Bayes QL dispersions.} 
\code{df.prior} is a vector of length \code{nrow(y)} if \code{robust=TRUE}, otherwise it has length 1.
\code{var.prior} is a vector of length \code{nrow(y)} if \code{abundance.trend=TRUE}, otherwise it has length 1.

\code{glmQFTest} produce an object of class \code{DGELRT} with the same components as produced by \code{\link{glmLRT}}, except that the \code{table$LR} column becomes \code{table$F} and contains quasi-likelihood F-statistics.
It also stores \code{df.total}, a numeric vector containing the denominator degrees of freedom for the F-test, equal to \code{df.prior + df.residual.zeros}.
}

\references{
Chen Y, Lun ATL, and Smyth, GK (2016).
From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline.
\emph{F1000Research} 5, 1438.
\url{http://f1000research.com/articles/5-1438}

Lun, ATL, Chen, Y, and Smyth, GK (2016).
It's DE-licious: a recipe for differential expression analyses of RNA-seq experiments using quasi-likelihood methods in edgeR.
\emph{Methods in Molecular Biology} 1418, 391-416.
\url{http://www.statsci.org/smyth/pubs/QLedgeRPreprint.pdf} (Preprint 8 April 2015)

Lund, SP, Nettleton, D, McCarthy, DJ, and Smyth, GK (2012).
Detecting differential expression in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates.
\emph{Statistical Applications in Genetics and Molecular Biology} Volume 11, Issue 5, Article 8.
\url{http://www.statsci.org/smyth/pubs/QuasiSeqPreprint.pdf}

Lun, ATL, and Smyth, GK (2017).
No counts, no variance: allowing for loss of degrees of freedom when assessing biological variability from RNA-seq data.
\emph{Statistical Applications in Genetics and Molecular Biology} 16(2), 83-93.
\url{https://doi.org/10.1515/sagmb-2017-0010}

Phipson, B, Lee, S, Majewski, IJ, Alexander, WS, and Smyth, GK (2016).
Robust hyperparameter estimation protects against hypervariable genes and improves power to detect differential expression.
\emph{Annals of Applied Statistics} 10, 946-963.
\url{https://projecteuclid.org/euclid.aoas/1469199900}
}

\author{Yunshun Chen, Aaron Lun, Davis McCarthy and Gordon Smyth}

\seealso{
\code{\link{topTags}} displays results from \code{glmQLFTest}.

\code{\link{plotQLDisp}} can be used to visualize the distribution of QL dispersions after EB shrinkage from \code{glmQLFit}.

The \code{QuasiSeq} package gives an alternative implementation of the Lund et al (2012) methods.
}

\examples{
nlibs <- 4
ngenes <- 1000
dispersion.true <- 1/rchisq(ngenes, df=10)
design <- model.matrix(~factor(c(1,1,2,2)))

# Generate count data
y <- rnbinom(ngenes*nlibs,mu=20,size=1/dispersion.true)
y <- matrix(y,ngenes,nlibs)
d <- DGEList(y)
d <- calcNormFactors(d)

# Fit the NB GLMs with QL methods
d <- estimateDisp(d, design)
fit <- glmQLFit(d, design)
results <- glmQLFTest(fit)
topTags(results)
fit <- glmQLFit(d, design, robust=TRUE)
results <- glmQLFTest(fit)
topTags(results)
fit <- glmQLFit(d, design, abundance.trend=FALSE)
results <- glmQLFTest(fit)
topTags(results)
}

\concept{Model fit}
\concept{Dispersion estimation}
\concept{Differential expression}

