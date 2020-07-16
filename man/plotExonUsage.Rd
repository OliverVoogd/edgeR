\name{plotExonUsage}
\alias{plotExonUsage}

\title{Create a Plot of Exon Usage from Exon-Level Count Data}

\description{Create a plot of exon usage for a given gene by plotting the (un)transformed counts for each exon, coloured by experimental group.}

\usage{
plotExonUsage(y, geneID, group=NULL, transform="none", counts.per.million=TRUE,
              legend.coords=NULL, \dots)
}
\arguments{ 

\item{y}{either a matrix of exon-level counts, a list containing a matrix of counts for each exon or a \code{DGEList} object with (at least) elements \code{counts} (table of counts summarized at the exon level) and \code{samples} (data frame containing information about experimental group, library size and normalization factor for the library size). Each row of \code{y} should represent one exon.}

\item{geneID}{character string giving the name of the gene for which exon usage is to be plotted.} 

\item{group}{factor supplying the experimental group/condition to which each sample (column of \code{y}) belongs. If \code{NULL} (default) the function will try to extract if from \code{y}, which only works if \code{y} is a \code{DGEList} object.}

\item{transform}{character, supplying the method of transformation to be applied to the exon counts, if any. Options are \code{"none"} (original counts are preserved), \code{"sqrt"} (square-root transformation) and \code{"log2"} (log2 transformation). Default is \code{"none"}.}

\item{counts.per.million}{logical, if \code{TRUE} then counts per million (as determined from total library sizes) will be plotted for each exon, if \code{FALSE} the raw read counts will be plotted. Using counts per million effectively normalizes for different read depth among the different samples, which can make the exon usage plots easier to interpret.}

\item{legend.coords}{optional vector of length 2 giving the x- and y-coordinates of the legend on the plot. If \code{NULL} (default), the legend will be automatically placed near the top right corner of the plot.}

\item{\dots}{optional further arguments to be passed on to \code{plot}.}
}

\value{\code{plotExonUsage} (invisibly) returns the transformed matrix of counts for the gene being plotted and produces a plot to the current device.}

\details{
This function produces a simple plot for comparing exon usage between different experimental conditions for a given gene.
}

\author{Davis McCarthy, Gordon Smyth}

\examples{
# generate exon counts from NB, create list object
y<-matrix(rnbinom(40,size=1,mu=10),nrow=10)
rownames(y) <- rep(c("gene.1","gene.2"), each=5)
d<-DGEList(counts=y,group=rep(1:2,each=2))
plotExonUsage(d, "gene.1")
}

\seealso{
\code{\link{spliceVariants}} for methods to detect genes with evidence for alternative exon usage.
}

\concept{Differential exon usage}
