\title{expandAsMatrix}
\name{expandAsMatrix}
\alias{expandAsMatrix}
\description{
Expand scalar or vector to a matrix.
}
\usage{
expandAsMatrix(x, dim=NULL, byrow=TRUE)
}
\arguments{
  \item{x}{scalar, vector, matrix or CompressedMatrix.}
  \item{dim}{integer vector of length 2 specifying the required dimensions of the output matrix.}
  \item{byrow}{logical scalar specifying if matrix should be filled by columns or rows for a vector \code{x}.}
  
}

\details{
This function expands a scalar, row/column vector or CompressedMatrix to be a matrix of dimensions \code{dim}.
It is used internally in edgeR to convert offsets, weights and other values to a matrix for consistent handling.
If \code{dim} is \code{NULL}, the function is equivalent to calling \code{as.matrix(x)}.

If \code{x} is a vector, its length must match one of the output dimensions.
The matrix will then be filled by repeating the matrix across the matching dimension.
For example, if \code{length(x)==dim[1]}, the matrix will be filled such that each row contains \code{x}.
If both dimensions match, filling is determined by \code{byrow}, with filling by rows as the default.

If \code{x} CompressedMatrix object, the size of any non-repeated dimensions must be consistent with corresponding output dimension.
The \code{byrow} argument will be ignored as the repeat specifications will dictate how expansion should be performed.
See \code{?\link{CompressedMatrix}} for more details.
}

\value{
Numeric matrix of dimension \code{dim}.
}
\author{Gordon Smyth}
\examples{
expandAsMatrix(1:3,c(4,3))
expandAsMatrix(1:4,c(4,3))
}
