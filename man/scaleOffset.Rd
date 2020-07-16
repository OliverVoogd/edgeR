\name{scaleOffset}
\alias{scaleOffset}
\alias{scaleOffset.DGEList}
\alias{scaleOffset.default}

\title{Scale offsets}

\description{Ensures scale of offsets are consistent with library sizes.}

\usage{
\method{scaleOffset}{DGEList}(y, offset, \dots)
\method{scaleOffset}{default}(y, offset, \dots)
}
\arguments{ 
\item{y}{numeric vector or matrix of counts, or a \code{DGEList} object.}
\item{offset}{numeric vector or matrix of offsets to be scaled. If a vector, its length must equal to the length of \code{y} or the number of columns of \code{y}. If a matrix, its dimension must equal to the dimension of \code{y}.}
\item{\dots}{other arguments that are not currently used.}
}

\value{
\code{scaleOffset.default} returns a numeric vector if \code{offset} is a vector or a matrix if \code{offset} is a matrix.
\code{scaleOffset.DGEList} computes the scaled offests and store them in the \code{offset} component of the input \code{DGEList} object.
}

\details{
\code{scaleOffset} ensures that the scale of offsets are consistent with library sizes.
This is done by ensuring that the mean offset for each gene is the same as the mean log-library size.
The length or dimensions of \code{offset} should be consistent with the number of libraries in \code{y}.
}

\author{Aaron Lun, Yunshun Chen}

\examples{
y <- matrix(rnbinom(40,size=1,mu=100),10,4)
offset <- rnorm(4)
scaleOffset(y, offset)
}

\concept{Normalization}

