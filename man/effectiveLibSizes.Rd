\name{effectiveLibSizes}
\alias{effectiveLibSizes}
\alias{effectiveLibSizes.default}
\alias{effectiveLibSizes.DGEList}
\alias{effectiveLibSizes.DGEGLM}
\alias{effectiveLibSizes.DGELRT}

\title{Effective Library Sizes}

\description{Extract effective (normalized) library sizes.}

\usage{
\method{effectiveLibSizes}{default}(y, log = FALSE, \dots)
}

\arguments{ 
\item{y}{
  a object of class \code{DGEList}, \code{DGEGLM} or \code{DGELRT}.
  Alternatively a numeric matrix or an object that can be coerced to a numeric matrix.
}
\item{log}{
  logical, if \code{TRUE} then the library sizes are return on the natural log scale.
}
\item{\dots}{
  other arguments are not currently used.
}
}

\value{
A numeric matrix of effective (normalized) library sizes.
If \code{log=TRUE}, then natural log values are returned, equal to library size offsets for a NB log-linear model.
}

\details{
This function extracts normalized library sizes, equal to the original library sizes multiplied by the corresponding normalization factors, from an edgeR data object or fitted model object.

If the object contains a row-specific offsets (i.e., a non-sparse matrix of offsets), then the offsets for the first row are returned.
}

\author{Gordon Smyth}

\seealso{
\code{\link{calcNormFactors}}
}

\examples{
ngenes <- 100
nsamples <- 4
y <- DGEList(counts=matrix(rnbinom(ngenes*nsamples,size=1,mu=10),ngenes,nsamples))
y <- calcNormFactors(y)
data.frame(y$samples, eff.lib.size=effectiveLibSizes(y))
}

\concept{Normalization}
