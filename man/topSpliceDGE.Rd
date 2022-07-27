\title{Top table of differentially spliced genes or exons}
\name{topSpliceDGE}
\alias{topSpliceDGE}
\description{
Top table ranking the most differentially spliced genes or exons.
}
\usage{
topSpliceDGE(lrt, test="Simes", number=10, FDR=1)
}
\arguments{
  \item{lrt}{\code{DGELRT} object produced by \code{diffSpliceDGE}.}
  \item{test}{character string, possible values are \code{"Simes"}, \code{"gene"} or \code{"exon"}.
    \code{"exon"} gives exon-level tests for each exon.
    \code{"gene"} gives gene-level tests for each gene.
    \code{"Simes"} gives genewise p-values derived from the exon-level tests after Simes adjustment for each gene.}
  \item{number}{integer, maximum number of rows to output.}
  \item{FDR}{numeric, only show exons or genes with false discovery rate less than this cutoff.}
}

\details{
Ranks genes or exons by evidence for differential splicing.
The exon-level tests test for differences between each exon and all the exons for the same gene.
The gene-level tests test for any differences in exon usage between experimental conditions.

The Simes method processes the exon-level p-values to give an overall call of differential splicing for each gene.
It returns the minimum Simes-adjusted p-values for each gene.

The gene-level tests are likely to be powerful for genes in which several exons are differentially splices.
The Simes p-values is likely to be more powerful when only a minority of the exons for a gene are differentially spliced.
The exon-level tests are not recommended for formal error rate control.
}

\value{A data.frame with any annotation columns found in \code{lrt} plus the following columns
  \item{NExons}{number of exons if \code{test="Simes"} or \code{"gene"}}
  \item{Gene.Exon}{exon annotation if \code{test="exon"}}
  \item{logFC}{log-fold change of one exon vs all the exons for the same gene (if \code{test="exon"})}
  \item{exon.LR}{LR-statistics for exons (if \code{test="exon"} and the object for \code{diffSpliceDGE} was produced by \code{glmFit})}
  \item{exon.F}{F-statistics for exons (if \code{test="exon"} and the object for \code{diffSpliceDGE} was produced by \code{glmQLFit})}
  \item{gene.LR}{LR-statistics for genes (if \code{test="gene"} and the object for \code{diffSpliceDGE} was produced by \code{glmFit})}
  \item{gene.F}{F-statistics for genes (if \code{test="gene"} and the object for \code{diffSpliceDGE} was produced by \code{glmQLFit})}
  \item{P.Value}{p-value}
  \item{FDR}{false discovery rate}
}

\author{Yunshun Chen and Gordon Smyth}

\seealso{ \code{\link{diffSpliceDGE}}.}
