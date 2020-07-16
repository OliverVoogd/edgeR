\name{DGELRT-class}
\docType{class}
\alias{DGELRT-class}
\alias{show,DGELRT-method}
\title{Digital Gene Expression Likelihood Ratio Test data and results - class}

\description{
A list-based S4 class for storing results of a GLM-based differential expression analysis for DGE data.
}

\section{List Components}{
For objects of this class, rows correspond to genomic features and columns to statistics associated with the differential expression analysis.
The genomic features are called genes, but in reality might correspond to transcripts, tags, exons etc.

Objects of this class contain the following list components:
\describe{
  \item{\code{table}:}{ data frame containing the log-concentration (i.e. expression level), the log-fold change in expression between the two groups/conditions and the exact p-value for differential expression, for each gene.}
  \item{\code{coefficients.full}:}{ matrix containing the coefficients
  computed from fitting the full model (fit using \code{glmFit} and a
  given design matrix) to each gene in the dataset.}
  \item{\code{coefficients.null}:}{ matrix containing the coefficients
  computed from fitting the null model to each gene in the
  dataset.  The null model is the model to which the full model is
  compared, and is fit using \code{glmFit} and dropping selected
  column(s) (i.e. coefficient(s)) from the design matrix for the full model.}
  \item{\code{design}:}{ design matrix for the full model from the likelihood
  ratio test.}
  \item{\code{\dots}:}{ if the argument \code{y} to \code{glmLRT} (which
  produces the \code{DGELRT} object) was itself a \code{DGEList} object, then
  the \code{DGELRT} will contain all of the elements of \code{y},
  except for the table of counts and the table of pseudocounts.}
}
}

\section{Methods}{
This class inherits directly from class \code{list}, so \code{DGELRT} objects can be manipulated as if they were ordinary lists.
However they can also be treated as if they were matrices for the purposes of subsetting.

The dimensions, row names and column names of a \code{DGELRT} object are defined by those of \code{table}, see \code{\link{dim.DGELRT}} or \code{\link{dimnames.DGELRT}}.

\code{DGELRT} objects can be subsetted, see \code{\link{subsetting}}.

\code{DGELRT} objects also have a \code{show} method so that printing produces a compact summary of their contents.
}

\author{edgeR team. First created by Davis McCarthy}

\seealso{
Other classes defined in edgeR are \code{\link{DGEList-class}}, \code{\link{DGEExact-class}}, \code{\link{DGEGLM-class}}, \code{\link{TopTags-class}}
}

\concept{edgeR classes}

