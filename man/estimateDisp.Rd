\name{estimateDisp}
\alias{estimateDisp}
\alias{estimateDisp.DGEList}
\alias{estimateDisp.SummarizedExperiment}
\alias{estimateDisp.default}

\title{Estimate Common, Trended and Tagwise Negative Binomial dispersions by weighted likelihood empirical Bayes}

\description{
Maximizes the negative binomial likelihood to give the estimate of the common, trended and tagwise dispersions across all tags.
}

\usage{
\S3method{estimateDisp}{DGEList}(y, design=NULL, prior.df=NULL, trend.method="locfit", tagwise=TRUE,
            span=NULL, min.row.sum=5, grid.length=21, grid.range=c(-10,10), robust=FALSE, 
            winsor.tail.p=c(0.05,0.1), tol=1e-06, ...)
\S3method{estimateDisp}{SummarizedExperiment}(y, design=NULL, prior.df=NULL, trend.method="locfit", tagwise=TRUE,
            span=NULL, min.row.sum=5, grid.length=21, grid.range=c(-10,10), robust=FALSE, 
            winsor.tail.p=c(0.05,0.1), tol=1e-06, ...)
\S3method{estimateDisp}{default}(y, design=NULL, group=NULL, lib.size=NULL, offset=NULL, prior.df=NULL,
            trend.method="locfit", tagwise=TRUE, span=NULL, min.row.sum=5, grid.length=21, 
            grid.range=c(-10,10), robust=FALSE, winsor.tail.p=c(0.05,0.1), tol=1e-06, weights=NULL, ...)
}

\arguments{
\item{y}{matrix of counts, or a \code{DGEList} object, or a \code{SummarizedExperiment} object.}
\item{design}{numeric design matrix. Defaults to \code{model.matrix(~group)} if \code{group} is specified and otherwise to a single column of ones.}
\item{prior.df}{prior degrees of freedom. It is used in calculating \code{prior.n}.}
\item{trend.method}{method for estimating dispersion trend. Possible values are \code{"locfit"} (default), \code{"none"}, \code{"movingave"}, \code{"loess"} and \code{"locfit.mixed"}, which uses a polynomial of degree 1 for lowly expressed genes.}
\item{tagwise}{logical, should the tagwise dispersions be estimated?}
\item{span}{width of the smoothing window, as a proportion of the data set.}
\item{min.row.sum}{numeric scalar giving a value for the filtering out of low abundance tags. Only tags with total sum of counts above this value are used. Low abundance tags can adversely affect the dispersion estimation, so this argument allows the user to select an appropriate filter threshold for the tag abundance.}
\item{grid.length}{the number of points on which the interpolation is applied for each tag.}
\item{grid.range}{the range of the grid points around the trend on a log2 scale.}
\item{robust}{logical, should the estimation of \code{prior.df} be robustified against outliers?}
\item{winsor.tail.p}{numeric vector of length 1 or 2, giving left and right tail proportions of the deviances to Winsorize when estimating \code{prior.df}.}
\item{tol}{the desired accuracy, passed to \code{\link{optimize}}}
\item{group}{vector or factor giving the experimental group/condition for each library. Defaults to a vector of ones with length equal to the number of libraries.}
\item{lib.size}{numeric vector giving the total count (sequence depth) for each library.}
\item{offset}{offset matrix for the log-linear model, as for \code{\link{glmFit}}. Defaults to the log-effective library sizes.}
\item{weights}{optional numeric matrix giving observation weights}
\item{\dots}{other arguments that are not currently used.}
}

\value{
\code{estimateDisp.DGEList} adds the following components to the input \code{DGEList} object:
	\item{design}{the design matrix.}
	\item{common.dispersion}{estimate of the common dispersion.}
	\item{trended.dispersion}{estimates of the trended dispersions.}
	\item{tagwise.dispersion}{tagwise estimates of the dispersion parameter if \code{tagwise=TRUE}.}
	\item{AveLogCPM}{numeric vector giving log2(AveCPM) for each row of \code{y}.}
	\item{trend.method}{method for estimating dispersion trend as given in the input.}
	\item{prior.df}{prior degrees of freedom. It is a vector when robust method is used.}
	\item{prior.n}{estimate of the prior weight, i.e. the smoothing parameter that indicates the weight to put on the common likelihood compared to the individual tag's likelihood.}
	\item{span}{width of the smoothing window used in estimating dispersions.}

\code{estimateDisp.SummarizedExperiment} converts the input \code{SummarizedExperiment} object into a \code{DGEList} object, and then calls \code{estimateDisp.DGEList}. The output is a \code{DGEList} object.

\code{estimateDisp.default} returns a list containing \code{common.dispersion}, \code{trended.dispersion}, \code{tagwise.dispersion} (if \code{tagwise=TRUE}), \code{span}, \code{prior.df} and \code{prior.n}.
}

\details{
This function calculates a matrix of likelihoods for each tag at a set of dispersion grid points, and then applies weighted likelihood empirical Bayes method to obtain posterior dispersion estimates. If there is no design matrix, it calculates the quantile conditional likelihood for each tag and then maximizes it. In this case, it is similar to the function \code{estimateCommonDisp} and \code{estimateTagwiseDisp}. If a design matrix is given, it calculates the adjusted profile log-likelihood for each tag and then maximizes it. In this case, it is similar to the functions \code{estimateGLMCommonDisp}, \code{estimateGLMTrendedDisp} and \code{estimateGLMTagwiseDisp}.

Note that the terms `tag' and `gene' are synonymous here.
}

\note{
The \code{estimateDisp} function doesn't give exactly the same estimates as the traditional calling sequences.
}

\references{
Chen, Y, Lun, ATL, and Smyth, GK (2014).
Differential expression analysis of complex RNA-seq experiments using edgeR.
In: \emph{Statistical Analysis of Next Generation Sequence Data},
Somnath Datta and Daniel S. Nettleton (eds), Springer, New York, pages 51-74.
\url{http://www.statsci.org/smyth/pubs/edgeRChapterPreprint.pdf}

Phipson, B, Lee, S, Majewski, IJ, Alexander, WS, and Smyth, GK (2016).
Robust hyperparameter estimation protects against hypervariable genes and improves power to detect differential expression.
\emph{Annals of Applied Statistics} 10, 946-963.
\url{http://projecteuclid.org/euclid.aoas/1469199900}
}

\author{Yunshun Chen, Gordon Smyth}
\examples{
# True dispersion is 1/5=0.2
y <- matrix(rnbinom(1000, mu=10, size=5), ncol=4)
group <- factor(c(1,1,2,2))
design <- model.matrix(~group)
d <- DGEList(counts=y, group=group)
d1 <- estimateDisp(d)
d2 <- estimateDisp(d, design)
}

\seealso{
\code{\link{estimateCommonDisp}}, \code{\link{estimateTagwiseDisp}}, \code{\link{estimateGLMCommonDisp}}, \code{\link{estimateGLMTrendedDisp}}, \code{\link{estimateGLMTagwiseDisp}}
}

\concept{Dispersion estimation}
