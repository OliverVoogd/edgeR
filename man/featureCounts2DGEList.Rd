\name{featureCounts2DGEList}
\alias{featureCounts2DGEList}

\title{Convert featureCounts object to a DGEList}

\description{
Converts the list output by \code{Rsubread::featureCounts} to a DGEList object.
}

\usage{
featureCounts2DGEList(x)
}

\arguments{
  \item{x}{a list produced by \code{Rsubread::featureCounts.}}
}

\details{
Rsubread's \code{featureCounts} function counts reads by features (typically exons or genomic intervals) or meta-features (typically genes).
It may also count reads crossing exon-exon junctions or reads internal to exons.
This function assembles all the read counts output by \code{featureCounts} into an DGEList, ensuring unique row.names where appropriate.
The proportion of reads assigned to features is also stored.
}

\value{A \code{\link[edgeR:DGEList-class]{DGEList}} data object.}

\author{Gordon Smyth.}

\seealso{\code{\link[edgeR:DGEList-class]{DGEList-class}}}

\concept{edgeR classes}
