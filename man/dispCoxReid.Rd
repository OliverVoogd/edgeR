\name{dispCoxReid}
\alias{dispCoxReid}
\alias{dispDeviance}
\alias{dispPearson}

\title{Estimate Common Dispersion for Negative Binomial GLMs}

\description{
Estimate a common dispersion parameter across multiple negative binomial generalized linear models.
}

\usage{
dispCoxReid(y, design=NULL, offset=NULL, weights=NULL, AveLogCPM=NULL, interval=c(0,4),
            tol=1e-5, min.row.sum=5, subset=10000)
dispDeviance(y, design=NULL, offset=NULL, interval=c(0,4), tol=1e-5, min.row.sum=5,
            subset=10000, AveLogCPM=NULL, robust=FALSE, trace=FALSE)
dispPearson(y, design=NULL, offset=NULL, min.row.sum=5, subset=10000,
            AveLogCPM=NULL, tol=1e-6, trace=FALSE, initial.dispersion=0.1)
}

\arguments{
\item{y}{numeric matrix of counts. A glm is fitted to each row.}

\item{design}{numeric design matrix, as for \code{\link{glmFit}}.}

\item{offset}{numeric vector or matrix of offsets for the log-linear models, as for \code{\link{glmFit}}. Defaults to \code{log(colSums(y))}.}

\item{weights}{optional numeric matrix giving observation weights}

\item{AveLogCPM}{numeric vector giving average log2 counts per million.}

\item{interval}{numeric vector of length 2 giving minimum and maximum allowable values for the dispersion, passed to \code{optimize}.}

\item{tol}{the desired accuracy, see \code{optimize} or \code{uniroot}.}

\item{min.row.sum}{integer. Only rows with at least this number of counts are used.}

\item{subset}{integer, number of rows to use in the calculation.  Rows used are chosen evenly spaced by AveLogCPM.}

\item{trace}{logical, should iteration information be output?}

\item{robust}{logical, should a robust estimator be used?}

\item{initial.dispersion}{starting value for the dispersion}
}

\value{
Numeric vector of length one giving the estimated common dispersion.
}

\details{
These are low-level (non-object-orientated) functions called by \code{estimateGLMCommonDisp}.

\code{dispCoxReid} maximizes the Cox-Reid adjusted profile likelihood (Cox and Reid, 1987).
\code{dispPearson} sets the average Pearson goodness of fit statistics to its (asymptotic) expected value.
This is also known as the \emph{pseudo-likelihood} estimator.
\code{dispDeviance} sets the average residual deviance statistic to its (asymptotic) expected values.
This is also known as the \emph{quasi-likelihood} estimator.

Robinson and Smyth (2008) and McCarthy et al (2011) showed that the Pearson (pseudo-likelihood) estimator typically under-estimates the true dispersion.
It can be seriously biased when the number of libraries (\code{ncol(y)} is small.
On the other hand, the deviance (quasi-likelihood) estimator typically over-estimates the true dispersion when the number of libraries is small.
Robinson and Smyth (2008) and McCarthy et al (2011) showed the Cox-Reid estimator to be the least biased of the three options.

\code{dispCoxReid} uses \code{optimize} to maximize the adjusted profile likelihood.
\code{dispDeviance} uses \code{uniroot} to solve the estimating equation.
The robust options use an order statistic instead the mean statistic, and have the effect that a minority of genes with very large (outlier) dispersions should have limited influence on the estimated value.
\code{dispPearson} uses a globally convergent Newton iteration.
}

\references{
Cox, DR, and Reid, N (1987). Parameter orthogonality and approximate conditional inference.
\emph{Journal of the Royal Statistical Society Series B} 49, 1-39.

Robinson MD and Smyth GK (2008). Small-sample estimation of negative
binomial dispersion, with applications to SAGE data.
\emph{Biostatistics}, 9, 321-332

McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
\emph{Nucleic Acids Research}.
\url{http://nar.oxfordjournals.org/content/early/2012/02/06/nar.gks042}
(Published online 28 January 2012)
}

\author{Gordon Smyth}
\examples{
ngenes <- 100
nlibs <- 4
y <- matrix(rnbinom(ngenes*nlibs,mu=10,size=10),nrow=ngenes,ncol=nlibs)
group <- factor(c(1,1,2,2))
lib.size <- rowSums(y)
design <- model.matrix(~group)
disp <- dispCoxReid(y, design, offset=log(lib.size), subset=100)
}

\seealso{
\code{\link{estimateGLMCommonDisp}}, \code{\link{optimize}}, \code{\link{uniroot}}
}
