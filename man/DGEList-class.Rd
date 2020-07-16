\name{DGEList-class}
\alias{DGEList-class}
\docType{class}
\title{Digital Gene Expression data - class}

\description{
A list-based S4 class for storing read counts and associated information from digital gene expression or sequencing technologies.
}

\section{List Components}{
For objects of this class, rows correspond to genomic features and columns to samples.
The genomic features are called genes, but in reality might correspond to transcripts, tags, exons etc.
Objects of this class contain the following essential list components:
\describe{
  \item{\code{counts}:}{ numeric matrix of read counts, one row for each gene and one column for each sample.}
  \item{\code{samples}:}{ data.frame with a row for each sample and columns \code{group}, \code{lib.size} and \code{norm.factors} containing the group labels, library sizes and normalization factors.
  Other columns can be optionally added to give more detailed sample information.}
}
Optional components include:
\describe{
  \item{\code{genes}:}{ data.frame giving annotation information for each gene.  Same number of rows as \code{counts}.}
  \item{\code{AveLogCPM}:}{ numeric vector giving average log2 counts per million for each gene.}
  \item{\code{common.dispersion}:}{ numeric scalar giving the overall dispersion estimate.}
  \item{\code{trended.dispersion}:}{ numeric vector giving trended dispersion estimates for each gene.}
  \item{\code{tagwise.dispersion}:}{ numeric vector giving tagwise dispersion estimates for each gene (note that `tag' and `gene' are synonymous here).}
  \item{\code{offset}:}{ numeric matrix of same size as \code{counts} giving offsets for use in log-linear models.}
}
}

\seealso{
\code{\link{DGEList}} constructs DGEList objects.
Other classes defined in edgeR are \code{\link{DGEExact-class}}, \code{\link{DGEGLM-class}}, \code{\link{DGELRT-class}}, \code{\link{TopTags-class}}
}

\section{Methods}{
This class inherits directly from class \code{list}, so \code{DGEList} objects can be manipulated as if they were ordinary lists.
However they can also be treated as if they were matrices for the purposes of subsetting.

The dimensions, row names and column names of a \code{DGEList} object are defined by those of \code{counts}, see \code{\link{dim.DGEList}} or \code{\link{dimnames.DGEList}}.

\code{DGEList} objects can be subsetted, see \code{\link[edgeR:subsetting]{subsetting}}.

\code{DGEList} objects also have a \code{show} method so that printing produces a compact summary of their contents.
}

\author{edgeR team. First created by Mark Robinson.}

\concept{edgeR classes}

