\name{estimateGLMCommonDisp}
\alias{estimateGLMCommonDisp}
\alias{estimateGLMCommonDisp.DGEList}
\alias{estimateGLMCommonDisp.default}

\title{Estimate Common Dispersion for Negative Binomial GLMs}

\description{
Estimates a common negative binomial dispersion parameter for a DGE dataset with a general experimental design.
}

\usage{
\S3method{estimateGLMCommonDisp}{DGEList}(y, design=NULL, method="CoxReid",
                      subset=10000, verbose=FALSE, \dots)
\S3method{estimateGLMCommonDisp}{default}(y, design=NULL, offset=NULL,
                      method="CoxReid", subset=10000, AveLogCPM=NULL,
                      verbose=FALSE, weights=NULL,\dots)
}

\arguments{
\item{y}{object containing read counts, as for \code{\link{glmFit}}.}

\item{design}{numeric design matrix, as for \code{\link{glmFit}}.}

\item{offset}{numeric vector or matrix of offsets for the log-linear models, as for \code{\link{glmFit}}.}

\item{method}{method for estimating the dispersion.
Possible values are \code{"CoxReid"}, \code{"Pearson"} or \code{"deviance"}.}

\item{subset}{maximum number of rows of \code{y} to use in the calculation.  Rows used are chosen evenly spaced by AveLogCPM using \code{\link{systematicSubset}}.}

\item{AveLogCPM}{numeric vector giving average log2 counts per million for each gene.}

\item{verbose}{logical, if \code{TRUE} estimated dispersion and BCV will be printed to standard output.}
\item{weights}{optional numeric matrix giving observation weights}

\item{\ldots}{other arguments are passed to lower-level functions.
See \code{\link{dispCoxReid}}, \code{\link{dispPearson}} and \code{\link{dispDeviance}} for details.}
}

\value{
The default method returns a numeric vector of length 1 containing the estimated common dispersion.

The \code{DGEList} method returns the same \code{DGEList} \code{y} as input but with \code{common.dispersion} as an added component.
The output object will also contain a component \code{AveLogCPM} if it was not already present in \code{y}.
}

\details{
This function calls \code{dispCoxReid}, \code{dispPearson} or \code{dispDeviance} depending on the \code{method} specified.
See \code{\link{dispCoxReid}} for details of the three methods and a discussion of their relative performance.
}


\references{
McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
\emph{Nucleic Acids Research} 40, 4288-4297.
\doi{10.1093/nar/gks042}
}

\author{Gordon Smyth, Davis McCarthy, Yunshun Chen}
\examples{
#  True dispersion is 1/size=0.1
y <- matrix(rnbinom(1000,mu=10,size=10),ncol=4)
d <- DGEList(counts=y,group=c(1,1,2,2))
design <- model.matrix(~group, data=d$samples)
d1 <- estimateGLMCommonDisp(d, design, verbose=TRUE)

#  Compare with classic CML estimator:
d2 <- estimateCommonDisp(d, verbose=TRUE)

#  See example(glmFit) for a different example
}

\seealso{
\code{\link{dispCoxReid}}, \code{\link{dispPearson}}, \code{\link{dispDeviance}}

\code{\link{estimateGLMTrendedDisp}} for trended dispersions or \code{\link{estimateGLMTagwiseDisp}} for genewise dispersions in the context of a generalized linear model.

\code{\link{estimateCommonDisp}} for the common dispersion or \code{\link{estimateTagwiseDisp}} for genewise dispersions in the context of a multiple group experiment (one-way layout).
}

\concept{Dispersion estimation}
