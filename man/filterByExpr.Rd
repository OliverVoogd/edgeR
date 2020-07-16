\name{filterByExpr}
\alias{filterByExpr}
\alias{filterByExpr.DGEList}
\alias{filterByExpr.SummarizedExperiment}
\alias{filterByExpr.default}

\title{Filter Genes By Expression Level}

\description{Determine which genes have sufficiently large counts to be retained in a statistical analysis.}

\usage{
\method{filterByExpr}{DGEList}(y, design = NULL, group = NULL, lib.size = NULL, \dots)
\method{filterByExpr}{SummarizedExperiment}(y, design = NULL, group = NULL, lib.size = NULL, \dots)
\method{filterByExpr}{default}(y, design = NULL, group = NULL, lib.size = NULL,
             min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7, \dots)
}

\arguments{ 
\item{y}{matrix of counts, or a \code{DGEList} object, or a \code{SummarizedExperiment} object.}
\item{design}{design matrix. Ignored if \code{group} is not \code{NULL}.}
\item{group}{vector or factor giving group membership for a oneway layout, if appropriate.}
\item{lib.size}{library size, defaults to \code{colSums(y)}.}
\item{min.count}{numeric. Minimum count required for at least some samples.}
\item{min.total.count}{numeric. Minimum total count required.}
\item{large.n}{integer. Number of samples per group that is considered to be \dQuote{large}.}
\item{min.prop}{numeric. Minimum proportion of samples in the smallest group that express the gene.}
\item{\dots}{any other arguments.
For the \code{DGEList} and \code{SummarizedExperiment} methods, other arguments will be passed to the default method.
For the default method, other arguments are not currently used.}
}

\details{
This function implements the filtering strategy that was intuitively described by Chen et al (2016).
Roughly speaking, the strategy keeps genes that have at least \code{min.count} reads in a worthwhile number samples.
More precisely, the filtering keeps genes that have count-per-million (CPM) above \emph{k} in \emph{n} samples, where \emph{k} is determined by \code{min.count} and by the sample library sizes and \emph{n} is determined by the design matrix.

\emph{n} is essentially the smallest group sample size or, more generally, the minimum inverse leverage of any fitted value.
If all the group sizes are larger than \code{large.n}, then this is relaxed slightly, but with \emph{n} always greater than \code{min.prop} of the smallest group size (70\% by default).

In addition, each kept gene is required to have at least \code{min.total.count} reads across all the samples.
}

\value{
Logical vector of length \code{nrow(y)} indicating which rows of \code{y} to keep in the analysis.
}

\author{Gordon Smyth}

\references{
Chen Y, Lun ATL, and Smyth, GK (2016).
From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline.
\emph{F1000Research} 5, 1438.
\url{http://f1000research.com/articles/5-1438}
}

\examples{\dontrun{
keep <- filterByExpr(y, design)
y <- y[keep,]
}}
