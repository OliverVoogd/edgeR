\name{estimateGLMTrendedDisp}
\alias{estimateGLMTrendedDisp}
\alias{estimateGLMTrendedDisp.DGEList}
\alias{estimateGLMTrendedDisp.default}

\title{Estimate Trended Dispersion for Negative Binomial GLMs}

\description{
Estimates the abundance-dispersion trend by Cox-Reid approximate profile likelihood.
}

\usage{
\S3method{estimateGLMTrendedDisp}{DGEList}(y, design=NULL, method="auto", \dots)
\S3method{estimateGLMTrendedDisp}{default}(y, design=NULL, offset=NULL, AveLogCPM=NULL,
                       method="auto", weights=NULL, \dots)
}

\arguments{ 
\item{y}{a matrix of counts or a \code{DGEList} object.)}
\item{design}{numeric design matrix, as for \code{\link{glmFit}}.}
\item{method}{method (low-level function) used to estimated the trended dispersions.
Possible values are \code{"auto"} (default, switch to \code{"bin.spline"} method if the number of genes is great than 200 and \code{"power"} method otherwise),\code{"bin.spline"}, \code{"bin.loess"} (which both result in a call to \code{dispBinTrend}), \code{"power"} (call to \code{dispCoxReidPowerTrend}), or \code{"spline"} (call to \code{dispCoxReidSplineTrend}).}
\item{offset}{numeric scalar, vector or matrix giving the linear model offsets, as for \code{\link{glmFit}}.}
\item{AveLogCPM}{numeric vector giving average log2 counts per million for each gene.}
\item{weights}{optional numeric matrix giving observation weights}
\item{\ldots}{other arguments are passed to lower-level functions \code{\link{dispBinTrend}}, \code{\link{dispCoxReidPowerTrend}} or \code{\link{dispCoxReidSplineTrend}}.}
}

\value{
When the input object is a \code{DGEList}, \code{estimateGLMTrendedDisp} produces a \code{DGEList} object, which contains the estimates of the trended dispersion parameter for the negative binomial model according to the method applied.

When the input object is a numeric matrix, it returns a vector of trended dispersion estimates calculated by one of the lower-level functions \code{dispBinTrend}, \code{dispCoxReidPowerTrend} and \code{dispCoxReidSplineTrend}.
}

\details{
Estimates the dispersion parameter for each gene with a trend that depends on the overall level of expression for that gene.
This is done for a DGE dataset for general experimental designs by using Cox-Reid approximate conditional inference for a negative binomial generalized linear model for each gene with the unadjusted counts and design matrix provided.

The function provides an object-orientated interface to lower-level functions.
}


\references{
Cox, DR, and Reid, N (1987). Parameter orthogonality and approximate conditional inference. \emph{Journal of the Royal Statistical Society Series B} 49, 1-39.

McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
\emph{Nucleic Acids Research} 40, 4288-4297.
\doi{10.1093/nar/gks042}
}

\author{Gordon Smyth, Davis McCarthy, Yunshun Chen}
\examples{
ngenes <- 250
nlibs <- 4
y <- matrix(rnbinom(ngenes*nlibs,mu=10,size=10),ngenes,nlibs)
d <- DGEList(counts=y,group=c(1,1,2,2),lib.size=c(1000:1003))
design <- model.matrix(~group, data=d$samples)
disp <- estimateGLMTrendedDisp(d, design, min.n=25, df=3)
plotBCV(disp)
}

\seealso{
\code{\link{dispBinTrend}}, \code{\link{dispCoxReidPowerTrend}} and \code{\link{dispCoxReidSplineTrend}} for details on how the calculations are done.
}

\concept{Dispersion estimation}
