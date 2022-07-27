\name{estimateTrendedDisp}
\alias{estimateTrendedDisp}
\alias{estimateTrendedDisp.DGEList}
\alias{estimateTrendedDisp.default}

\title{Estimate Empirical Bayes Trended Dispersion Values}

\description{
Estimates trended dispersion values by an empirical Bayes method.
}

\usage{
\S3method{estimateTrendedDisp}{DGEList}(y, method="bin.spline", df=5, span=2/3, ...)
\S3method{estimateTrendedDisp}{default}(y, group=NULL, lib.size=NULL, AveLogCPM=NULL, 
            method="bin.spline", df=5, span=2/3, ...)
}

\arguments{ 
\item{y}{matrix of counts or a \code{DGEList} object.}
\item{method}{method used to estimated the trended dispersions. Possible values are \code{"bin.spline"}, and \code{"bin.loess"}.}
\item{df}{integer giving the degrees of freedom of the spline function if \code{"bin.spline"} method is used, see \code{ns} in the splines package. Default is 5.}
\item{span}{scalar, passed to \code{loess} to determine the amount of smoothing for the loess fit when \code{"loess"} method is used. Default is \code{2/3}.}
\item{group}{vector or factor giving the experimental group/condition for each library.}
\item{lib.size}{numeric vector giving the total count (sequence depth) for each library.}
\item{AveLogCPM}{numeric vector giving average log2 counts per million for each tag}
\item{\dots}{other arguments that are not currently used.}
}

\details{
This function takes the binned common dispersion and abundance, and fits a smooth curve through these binned values using either natural cubic splines or loess. From this smooth curve it predicts the dispersion value for each gene based on the gene's overall abundance. This results in estimates for the NB dispersion parameter which have a dependence on the overall expression level of the gene, and thus have an abundance-dependent trend.
}

\value{
An object of class \code{DGEList} with the same components as for \code{\link{estimateCommonDisp}} plus the trended dispersion estimates for each gene.
}

\author{Yunshun Chen and Gordon Smyth}
\examples{
ngenes <- 1000
nlib <- 4
log2cpm <- seq(from=0,to=16,length=ngenes)
lib.size <- 1e7
mu <- 2^log2cpm * lib.size * 1e-6
dispersion <- 1/sqrt(mu) + 0.1
counts <- rnbinom(ngenes*nlib, mu=mu, size=1/dispersion)
counts <- matrix(counts,ngenes,nlib)
y <- DGEList(counts,lib.size=rep(lib.size,nlib))
y <- estimateCommonDisp(y)
y <- estimateTrendedDisp(y)
}

\seealso{
\code{\link{estimateCommonDisp}} estimates a common value for the dispersion parameter for all genes - should generally be run before \code{estimateTrendedDisp}.
}

\concept{Dispersion estimation}
