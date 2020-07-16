\title{Differential splicing plot}
\name{plotSpliceDGE}
\alias{plotSpliceDGE}
\description{
Plot relative log-fold changes by exons for the specified gene and highlight the significantly spliced exons.
}
\usage{
plotSpliceDGE(lrt, geneid=NULL, genecolname=NULL, rank=1L, FDR=0.05)
}
\arguments{
  \item{lrt}{\code{DGELRT} object produced by \code{diffSpliceDGE}.}
  \item{geneid}{character string, ID of the gene to plot.}
  \item{genecolname}{column name of \code{lrt$genes} containing gene IDs. Defaults to \code{lrt$genecolname}.}
  \item{rank}{integer, if \code{geneid=NULL} then this ranked gene will be plotted.}
  \item{FDR}{numeric, mark exons with false discovery rate less than this cutoff.}
}

\details{
Plot relative log2-fold-changes by exon for the specified gene.
The relative logFC is the difference between the exon's logFC and the overall logFC for the gene, as computed by \code{diffSpliceDGE}.
The significantly spliced individual exons are highlighted as red dots. 
The size of the red dots are weighted by its significance.
}

\value{A plot is created on the current graphics device.}
\author{Yunshun Chen, Yifang Hu and Gordon Smyth}
\seealso{
\code{\link{diffSpliceDGE}}, \code{\link{topSpliceDGE}}.
}

\concept{Differential exon usage}
