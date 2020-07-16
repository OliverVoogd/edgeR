\name{estimateCommonDisp}
\alias{estimateCommonDisp}
\alias{estimateCommonDisp.DGEList}
\alias{estimateCommonDisp.default}


\title{Estimate Common Negative Binomial Dispersion by Conditional Maximum Likelihood}

\description{
Maximizes the negative binomial conditional common likelihood to estimate a common dispersion value across all genes.
}

\usage{
\method{estimateCommonDisp}{DGEList}(y, tol=1e-06, rowsum.filter=5, verbose=FALSE, ...)
\method{estimateCommonDisp}{default}(y, group=NULL, lib.size=NULL, tol=1e-06, 
          rowsum.filter=5, verbose=FALSE, ...)
}

\arguments{
\item{y}{matrix of counts or a \code{DGEList} object.}
\item{tol}{the desired accuracy, passed to \code{\link{optimize}}.}
\item{rowsum.filter}{genes with total count (across all samples) below this value will be filtered out before estimating the dispersion.}
\item{verbose}{logical, if \code{TRUE} then the estimated dispersion and BCV will be printed to standard output.}
\item{group}{vector or factor giving the experimental group/condition for each library.}
\item{lib.size}{numeric vector giving the total count (sequence depth) for each library.}
\item{\dots}{other arguments that are not currently used.}
}

\value{
\code{estimateCommonDisp.DGEList} adds the following components to the input \code{DGEList} object:
	\item{common.dispersion}{estimate of the common dispersion.}
	\item{pseudo.counts}{numeric matrix of pseudo-counts.}
	\item{pseudo.lib.size}{the common library size to which the pseudo-counts have been adjusted.}
	\item{AveLogCPM}{numeric vector giving log2(AveCPM) for each row of \code{y}.}
\code{estimateCommonDisp.default} returns a numeric scalar of the common dispersion estimate.
}

\details{
Implements the conditional maximum likelihood (CML) method proposed by Robinson and Smyth (2008) for estimating a common dispersion parameter.
This method proves to be accurate and nearly unbiased even for small counts and small numbers of replicates.

The CML method involves computing a matrix of quantile-quantile normalized counts, called pseudo-counts.
The pseudo-counts are adjusted in such a way that the library sizes are equal for all samples, while preserving differences between groups and variability within each group.
The pseudo-counts are included in the output of the function, but are intended mainly for internal edgeR use.
}

\references{
Robinson MD and Smyth GK (2008).
Small-sample estimation of negative binomial dispersion, with applications to SAGE data.
\emph{Biostatistics}, 9, 321-332.
\url{http://biostatistics.oxfordjournals.org/content/9/2/321}
}

\author{Mark Robinson, Davis McCarthy, Gordon Smyth}
\examples{
# True dispersion is 1/5=0.2
y <- matrix(rnbinom(250*4,mu=20,size=5),nrow=250,ncol=4)
dge <- DGEList(counts=y,group=c(1,1,2,2))
dge <- estimateCommonDisp(dge, verbose=TRUE)
}

\seealso{
\code{\link{equalizeLibSizes}},
\code{\link{estimateTrendedDisp}},
\code{\link{estimateTagwiseDisp}}
}

\concept{Dispersion estimation}
