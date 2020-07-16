\name{cpm}
\alias{cpm}
\alias{cpm.DGEList}
\alias{cpm.SummarizedExperiment}
\alias{cpm.DGEGLM}
\alias{cpm.DGELRT}
\alias{cpm.default}
\alias{rpkm}
\alias{rpkm.DGEList}
\alias{rpkm.SummarizedExperiment}
\alias{rpkm.DGEGLM}
\alias{rpkm.DGELRT}
\alias{rpkm.default}
\alias{cpmByGroup}
\alias{cpmByGroup.DGEList}
\alias{cpmByGroup.SummarizedExperiment}
\alias{cpmByGroup.default}
\alias{rpkmByGroup}
\alias{rpkmByGroup.DGEList}
\alias{rpkmByGroup.SummarizedExperiment}
\alias{rpkmByGroup.default}

\title{Counts per Million or Reads per Kilobase per Million}

\description{Compute counts per million (CPM) or reads per kilobase per million (RPKM).}

\usage{
\method{cpm}{DGEList}(y, normalized.lib.sizes = TRUE,
       log = FALSE, prior.count = 2, \dots)
\method{cpm}{SummarizedExperiment}(y, normalized.lib.sizes = TRUE,
       log = FALSE, prior.count = 2, \dots)
\method{cpm}{DGEGLM}(y, log = FALSE, shrunk = TRUE, \dots)
\method{cpm}{default}(y, lib.size = NULL,
       log = FALSE, prior.count = 2, \dots)
\method{rpkm}{DGEList}(y, gene.length = NULL, normalized.lib.sizes = TRUE,
       log = FALSE, prior.count = 2, \dots)
\method{rpkm}{SummarizedExperiment}(y, gene.length = NULL, normalized.lib.sizes = TRUE,
       log = FALSE, prior.count = 2, \dots)
\method{rpkm}{DGEGLM}(y, gene.length, log = FALSE, shrunk = TRUE, \dots)
\method{rpkm}{default}(y, gene.length, lib.size = NULL,
       log = FALSE, prior.count = 2, \dots)
\method{cpmByGroup}{DGEList}(y, group = NULL, dispersion = NULL, \ldots)
\method{cpmByGroup}{SummarizedExperiment}(y, group = NULL, dispersion = NULL, \ldots)
\method{cpmByGroup}{default}(y, group = NULL, dispersion = 0.05,
       offset = NULL, weights = NULL, log = FALSE, prior.count = 2, \ldots)
\method{rpkmByGroup}{DGEList}(y, group = NULL, gene.length = NULL, dispersion = NULL, \ldots)
\method{rpkmByGroup}{SummarizedExperiment}(y, group = NULL, gene.length = NULL, dispersion = NULL, \ldots)
\method{rpkmByGroup}{default}(y, group = NULL, gene.length, dispersion = 0.05,
       offset = NULL, weights = NULL, log = FALSE, prior.count = 2, \ldots)
}

\arguments{ 
\item{y}{
  a matrix-like object containing counts.
  Can be a numeric matrix, a \code{DGEList} object, a \code{SummarizedExperiment} object with a \code{"counts"} assay, or any object that can be coerced to a matrix by \code{as.matrix}.
  For \code{cpm} and \code{rpkm}, it can also be a \code{DGEGLM} or \code{DGELRT} object.
}
\item{normalized.lib.sizes}{logical, use normalized library sizes?}
\item{lib.size}{library size, defaults to \code{colSums(y)}.}
\item{log}{logical, if \code{TRUE} then \code{log2} values are returned.}
\item{prior.count}{average count to be added to each observation to avoid taking log of zero. Used only if \code{log=TRUE}.}
\item{shrunk}{
  logical, if \code{TRUE} then the usual coefficients from the fitted object will be used, if \code{FALSE} then the unshrunk coefficients will be used.
}
\item{gene.length}{vector of length \code{nrow(y)} giving gene length in bases, or the name of the column \code{y$genes} containing the gene lengths.}
\item{group}{factor giving group membership for columns of \code{y}. Defaults to \code{y$sample$group} for the \code{DGEList} method and to a single level factor for the default method.}
\item{dispersion}{numeric vector of negative binomial dispersions.}
\item{offset}{numeric matrix of same size as \code{y} giving offsets for the log-linear models.  Can be a scalar or a vector of length \code{ncol(y)}, in which case it is expanded out to a matrix.}
\item{weights}{numeric vector or matrix of non-negative quantitative weights.
Can be a vector of length equal to the number of libraries, or a matrix of the same size as \code{y}.}
\item{\dots}{other arguments are not used.}
}

\value{
A numeric matrix of CPM or RPKM values, on the log2 scale if \code{log=TRUE}.
\code{cpm} and \code{rpkm} produce matrices of the same size as \code{y}.
If \code{y} was a data object, then observed values are returned.
If \code{y} was a fitted model object, then fitted values are returned.

\code{cpmByGroup} and \code{rpkmByGroup} produce matrices with a column for each level of \code{group}.
}

\details{
CPM or RPKM values are useful descriptive measures for the expression level of a gene.
By default, the normalized library sizes are used in the computation for \code{DGEList} objects but simple column sums for matrices.

If log-values are computed, then a small count, given by \code{prior.count} but scaled to be proportional to the library size, is added to \code{y} to avoid taking the log of zero.

The \code{rpkm} methods for \code{DGEList}, \code{DGEGLM} or \code{DGELRT} objects will try to find the gene lengths in a column of \code{y$genes} called \code{Length} or \code{length}.
Failing that, it will look for any column name containing \code{"length"} in any capitalization.

The \code{cpm} and \code{rpkm} methods for \code{DGEGLM} and \code{DGELRT} fitted model objects return fitted CPM or RPKM values.
If \code{shrunk=TRUE}, then the CPM or RPKM values will reflect the \code{prior.count} input to the original linear model fit.
If \code{shrunk=FALSE}, then the CPM or RPKM values will be computed with \code{prior.count=0}.
Note that the latter could result in taking the log of near-zero values if \code{log=TRUE}.

\code{cpmByGroup} and \code{rpkmByGroup} compute group average values on the unlogged scale.
}

\note{
\code{aveLogCPM(y)}, \code{rowMeans(cpm(y,log=TRUE))} and \code{log2(rowMeans(cpm(y))} all give slightly different results.
}

\author{Davis McCarthy, Gordon Smyth, Yunshun Chen, Aaron Lun}

\seealso{
\code{\link{aveLogCPM}}
}

\examples{
y <- matrix(rnbinom(20,size=1,mu=10),5,4)
cpm(y)

d <- DGEList(counts=y, lib.size=1001:1004)
cpm(d)
cpm(d,log=TRUE)

d$genes <- data.frame(Length=c(1000,2000,500,1500,3000))
rpkm(d)

cpmByGroup(d, group=c(1,1,2,2))

rpkmByGroup(d, group=c(1,1,2,2))
}

\concept{Data exploration}
