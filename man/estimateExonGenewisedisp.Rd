\name{estimateExonGenewiseDisp}
\alias{estimateExonGenewiseDisp}

\title{Estimate Genewise Dispersions from Exon-Level Count Data}

\description{Estimate a dispersion value for each gene from exon-level count data by collapsing exons into the genes to which they belong.}

\usage{
estimateExonGenewiseDisp(y, geneID, group=NULL)
}
\arguments{ 

\item{y}{either a matrix of exon-level counts or a \code{DGEList} object with (at least) elements \code{counts} (table of counts summarized at the exon level) and \code{samples} (data frame containing information about experimental group, library size and normalization factor for the library size). Each row of \code{y} should represent one exon.}

\item{geneID}{vector of length equal to the number of rows of \code{y}, which provides the gene identifier for each exon in \code{y}. These identifiers are used to group the relevant exons into genes for the gene-level analysis of splice variation.}

\item{group}{factor supplying the experimental group/condition to which each sample (column of \code{y}) belongs. If \code{NULL} (default) the function will try to extract if from \code{y}, which only works if \code{y} is a \code{DGEList} object.}

}

\value{\code{estimateExonGenewiseDisp} returns a vector of genewise dispersion estimates, one for each unique \code{geneID}.}

\details{
This function can be used to compute genewise dispersion estimates (for an experiment with a one-way, or multiple group, layout) from exon-level count data. \code{estimateCommonDisp} and \code{estimateTagwiseDisp} are used to do the computation and estimation, and the default arguments for those functions are used.
}

\author{Davis McCarthy, Gordon Smyth}

\examples{
# generate exon counts from NB, create list object
y<-matrix(rnbinom(40,size=1,mu=10),nrow=10)
d<-DGEList(counts=y,group=rep(1:2,each=2))
genes <- rep(c("gene.1","gene.2"), each=5)
estimateExonGenewiseDisp(d, genes)
}

\seealso{
\code{\link{estimateCommonDisp}} and related functions for estimating the dispersion parameter for the negative binomial model.
}

\concept{Dispersion estimation}
