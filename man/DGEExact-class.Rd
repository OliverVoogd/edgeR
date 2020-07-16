\name{DGEExact-class}
\docType{class}
\alias{DGEExact-class}
\alias{show,DGEExact-method}
\title{differential expression of Digital Gene Expression data - class}

\description{
A list-based S4 class for for storing results of a differential expression analysis for DGE data.
}

\section{List Components}{
For objects of this class, rows correspond to genomic features and columns to statistics associated with the differential expression analysis.
The genomic features are called genes, but in reality might correspond to transcripts, tags, exons etc.

Objects of this class contain the following list components:
\describe{
  \item{\code{table}:}{ data frame containing columns for the log2-fold-change, \code{logFC}, the average log2-counts-per-million, \code{logCPM}, and the two-sided p-value \code{PValue}.}
  \item{\code{comparison}:}{ vector giving the two experimental groups/conditions being compared.}
  \item{\code{genes}:}{ a data frame containing information about each gene (can be \code{NULL}).}
}
}

\section{Methods}{
This class inherits directly from class \code{list}, so \code{DGEExact} objects can be manipulated as if they were ordinary lists.
However they can also be treated as if they were matrices for the purposes of subsetting.

The dimensions, row names and column names of a \code{DGEExact} object are defined by those of \code{table}, see \code{\link{dim.DGEExact}} or \code{\link{dimnames.DGEExact}}.

\code{DGEExact} objects can be subsetted, see \code{\link{subsetting}}.

\code{DGEExact} objects also have a \code{show} method so that printing produces a compact summary of their contents.
}

\author{edgeR team. First created by Mark Robinson and Davis McCarthy}

\seealso{
Other classes defined in edgeR are \code{\link{DGEList-class}}, \code{\link{DGEGLM-class}}, \code{\link{DGELRT-class}}, \code{\link{TopTags-class}}
}

\concept{edgeR classes}

