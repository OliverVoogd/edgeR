\name{estimateTagwiseDisp}
\alias{estimateTagwiseDisp}
\alias{estimateTagwiseDisp.DGEList}
\alias{estimateTagwiseDisp.default}

\title{Estimate Empirical Bayes Tagwise Dispersion Values}

\description{
Estimates tagwise dispersion values by an empirical Bayes method based on weighted conditional maximum likelihood.
}

\usage{
\method{estimateTagwiseDisp}{DGEList}(y, prior.df=10, trend="movingave", span=NULL, method="grid", 
           grid.length=11, grid.range=c(-6,6), tol=1e-06, verbose=FALSE, ...)
\method{estimateTagwiseDisp}{default}(y, group=NULL, lib.size=NULL, dispersion, AveLogCPM=NULL, 
           prior.df=10, trend="movingave", span=NULL, method="grid", grid.length=11, 
           grid.range=c(-6,6), tol=1e-06, verbose=FALSE, ...)
}

\arguments{
\item{y}{matrix of counts or a \code{DGEList} object.}
\item{prior.df}{prior degrees of freedom.}
\item{trend}{method for estimating dispersion trend. Possible values are \code{"movingave"} (default), \code{"loess"} and \code{"none"}.}
\item{span}{width of the smoothing window, as a proportion of the data set.}
\item{method}{method for maximizing the posterior likelihood.
Possible values are \code{"grid"} (default) for interpolation on grid points or \code{"optimize"} to call the function of the same name.}
\item{grid.length}{for \code{method="grid"}, the number of points on which the interpolation is applied for each tag.}
\item{grid.range}{for \code{method="grid"}, the range of the grid points around the trend on a log2 scale.}
\item{tol}{for \code{method="optimize"}, the tolerance for Newton-Rhapson iterations.}
\item{verbose}{logical, if \code{TRUE} then diagnostic ouput is produced during the estimation process.}
\item{group}{vector or factor giving the experimental group/condition for each library.}
\item{lib.size}{numeric vector giving the total count (sequence depth) for each library.}
\item{dispersion}{common dispersion estimate, used as an initial estimate for the tagwise estimates.}
\item{AveLogCPM}{numeric vector giving average log2 counts per million for each tag}
\item{\dots}{other arguments that are not currently used.}
}

\details{
This function implements the empirical Bayes strategy proposed by Robinson and Smyth (2007)
for estimating the tagwise negative binomial dispersions.
The experimental design is assumed to be a oneway layout with one or more experimental groups.
The empirical Bayes posterior is implemented as a conditional likelihood with tag-specific weights.

The prior values for the dispersions are determined by a global trend.
The individual tagwise dispersions are then squeezed towards this trend.
The prior degrees of freedom determines the weight given to the prior.
The larger the prior degrees of freedom, the more the tagwise dispersions are squeezed towards the global trend.
If the number of libraries is large, the prior becomes less important and the tagwise dispersion are determined more by the individual tagwise data.

If \code{trend="none"}, then the prior dispersion is just a constant, the common dispersion.
Otherwise, the trend is determined by a moving average (\code{trend="movingave"}) or loess smoother applied to the tagwise conditional log-likelihood.
\code{method="loess"} applies a loess curve of degree 0 as implemented in \code{\link{loessByCol}}.

\code{method="optimize"} is not recommended for routine use as it is very slow.
It is included for testing purposes.

Note that the terms `tag' and `gene' are synonymous here. The function is only named `Tagwise' for historical reasons.
}

\value{
\code{estimateTagwiseDisp.DGEList} adds the following components to the input \code{DGEList} object:
	\item{prior.df}{prior degrees of freedom.}
	\item{prior.n}{estimate of the prior weight.}
	\item{tagwise.dispersion}{numeric vector of the tagwise dispersion estimates.}
	\item{span}{width of the smoothing window, in terms of proportion of the data set.}
\code{estimateTagwiseDisp.default} returns a numeric vector of the tagwise dispersion estimates.
}

\references{
Robinson, MD, and Smyth, GK (2007). Moderated statistical tests for
assessing differences in tag abundance. \emph{Bioinformatics} 23, 2881-2887.
\doi{10.1093/bioinformatics/btm453}
}

\author{Mark Robinson, Davis McCarthy, Yunshun Chen and Gordon Smyth}
\examples{
# True dispersion is 1/5=0.2
y <- matrix(rnbinom(250*4,mu=20,size=5),nrow=250,ncol=4)
dge <- DGEList(counts=y,group=c(1,1,2,2))
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)
}

\seealso{
\code{\link{estimateCommonDisp}} is usually run before \code{estimateTagwiseDisp}.

\code{\link{movingAverageByCol}} and \code{\link{loessByCol}} implement the moving average or loess smoothers.
}

\concept{Dispersion estimation}
