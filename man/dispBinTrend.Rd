\name{dispBinTrend}
\alias{dispBinTrend}

\title{Estimate Dispersion Trend by Binning for NB GLMs}

\description{
Estimate the abundance-dispersion trend by computing the common dispersion for bins of genes of similar AveLogCPM and then fitting a smooth curve.
}

\usage{
dispBinTrend(y, design=NULL, offset=NULL, df = 5, span=0.3, min.n=400,
             method.bin="CoxReid", method.trend="spline", AveLogCPM=NULL,
             weights=NULL, \dots)
}

\arguments{
\item{y}{numeric matrix of counts}
\item{design}{numeric matrix giving the design matrix for the GLM that is to be fit.}
\item{offset}{numeric scalar, vector or matrix giving the offset (in addition to the log of the effective library size) that is to be included in the NB GLM for the genes. If a scalar, then this value will be used as an offset for all genes and libraries. If a vector, it should be have length equal to the number of libraries, and the same vector of offsets will be used for each gene. If a matrix, then each library for each gene can have a unique offset, if desired. In \code{adjustedProfileLik} the \code{offset} must be a matrix with the same dimension as the table of counts.}
\item{df}{degrees of freedom for spline curve.}
\item{span}{span used for loess curve.}
\item{min.n}{minimim number of genes in a bins.}
\item{method.bin}{method used to estimate the dispersion in each bin. Possible values are \code{"CoxReid"}, \code{"Pearson"} or \code{"deviance"}.}
\item{method.trend}{type of curve to smooth the bins. Possible values are \code{"spline"} for a natural cubic regression spline or \code{"loess"} for a linear lowess curve.}
\item{AveLogCPM}{numeric vector giving average log2 counts per million for each gene}
\item{weights}{optional numeric matrix giving observation weights}
\item{\dots}{other arguments are passed to \code{estimateGLMCommonDisp}}
}

\value{
list with the following components:
\item{AveLogCPM}{numeric vector containing the overall AveLogCPM for each gene}
\item{dispersion}{numeric vector giving the trended dispersion estimate for each gene}
\item{bin.AveLogCPM}{numeric vector of length equal to \code{nbins} giving the average (mean) AveLogCPM for each bin}
\item{bin.dispersion}{numeric vector of length equal to \code{nbins} giving the estimated common dispersion for each bin}
}

\details{
Estimate a dispersion parameter for each of many negative binomial generalized linear models by computing the common dispersion for genes sorted into bins based on overall AveLogCPM.
A regression natural cubic splines or a linear loess curve is used to smooth the trend and extrapolate a value to each gene.

If there are fewer than \code{min.n} rows of \code{y} with at least one positive count, then one bin is used.
The number of bins is limited to 1000.
}

\references{
McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
\emph{Nucleic Acids Research} 40, 4288-4297.
\url{https://doi.org/10.1093/nar/gks042}
}

\author{Davis McCarthy and Gordon Smyth}
\examples{
ngenes <- 1000
nlibs <- 4
means <- seq(5,10000,length.out=ngenes)
y <- matrix(rnbinom(ngenes*nlibs,mu=rep(means,nlibs),size=0.1*means),nrow=ngenes,ncol=nlibs)
keep <- rowSums(y) > 0
y <- y[keep,]
group <- factor(c(1,1,2,2))
design <- model.matrix(~group) # Define the design matrix for the full model
out <- dispBinTrend(y, design, min.n=100, span=0.3)
with(out, plot(AveLogCPM, sqrt(dispersion)))
}

\seealso{
\code{\link{estimateGLMTrendedDisp}}
}
