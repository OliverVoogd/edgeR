\title{makeCompressedMatrix}
\name{makeCompressedMatrix}
\alias{makeCompressedMatrix}
\alias{CompressedMatrix}
\alias{[.CompressedMatrix}
\alias{[<-.CompressedMatrix}
\alias{length.CompressedMatrix}
\alias{as.matrix.CompressedMatrix}
\alias{Ops.CompressedMatrix}
\alias{cbind.CompressedMatrix}
\alias{rbind.CompressedMatrix}
\alias{dim.CompressedMatrix}
\alias{CompressedMatrix-class}

\description{
Construct a CompressedMatrix object from a scalar, vector or matrix.
}

\usage{
makeCompressedMatrix(x, dims, byrow=TRUE)

\method{dim}{CompressedMatrix}(x)
\method{length}{CompressedMatrix}(x)

\method{[}{CompressedMatrix}(x, i, j, drop=TRUE) 
\method{[}{CompressedMatrix}(x, i, j) <- value

\method{Ops}{CompressedMatrix}(e1, e2)

\method{rbind}{CompressedMatrix}(...)
\method{cbind}{CompressedMatrix}(...)

\method{as.matrix}{CompressedMatrix}(x, ...)
}

\arguments{
  \item{x}{For \code{makeCompressedMatrix}, a scalar, vector, matrix or CompressedMatrix object.
        For the S3 methods, a CompressedMatrix object.}
  \item{dims}{an integer vector indicating the matrix dimensions, ignored if \code{x} is already a matrix}
  \item{byrow}{logical. If \code{x} is a vector, should it be repeated across rows (default) or across columns?}
  \item{i, j}{subset indices to apply to \code{x}, which behave the same as indices for \code{matrix} subsetting}
  \item{drop}{logical, indicating whether or not to drop dimensions when subsetting to a single row/column}
  \item{value}{an array-like object or vector to be used to replace values in \code{x}}
  \item{e1, e2}{a CompressedMatrix object}
  \item{...}{multiple CompressedMatrix objects for \code{rbind} and \code{cbind}.
  Otherwise additional arguments that are ignored in \code{as.matrix}.}
}

\details{
The CompressedMatrix is used throughout edgeR to save space in storing offsets and (to a lesser extent) weights.
This is because, for routine analyses, offsets are the same for all genes so it makes little sense to expand it to the full dimensions of the count matrix.
Most functions will accept a CompressedMatrix as input to \code{offset} or \code{weights} arguments.
}

\section{Class construction}{
The \code{makeCompressedMatrix} function creates a CompressedMatrix object from \code{x}. 
The CompressedMatrix class inherits from a matrix and holds two logical scalar attributes \code{repeat.row} and \code{repeat.col}.
Each attribute specifies whether the values are to be repeated across rows and/or across columns.
This avoids the need to store redundant values in a full-sized matrix of dimensions \code{dim}, as would be done with \code{\link{expandAsMatrix}}.

To illustrate, consider that rows usually correspond to genes while columns usually correspond to libraries.
If we have a vector of library sizes, this will hold one unique value per library that is the same for all genes.
Thus, we should use \code{byrow=TRUE}, which will construct a CompressedMatrix object storing one row containing this vector.
Here, \code{repeat.row=TRUE} and \code{repeat.col=FALSE}, indicating that the row is to be repeated for all genes.

On the other hand, we may have a vector of gene-specific values that is the same for all libraries (e.g., dispersions).
In this case, we should use \code{byrow=FALSE} to construct the CompressedMatrix object.
This will store one column with \code{repeat.row=FALSE} and \code{repeat.col=TRUE}, indicating that the column should be repeated across libraries.

In cases where \code{x} is a scalar, \code{byrow} is ignored and both \code{repeat.row} and \code{repeat.col} will be \code{TRUE} by default.
If \code{x} is a matrix, both attributes will be \code{FALSE}.
If \code{x} is a CompressedMatrix, it will be returned without modification.
}

\section{Class methods}{
Subsetting of a CompressedMatrix object depends on the values of \code{repeat.row} and \code{repeat.col}.
If the rows are repeated, any subsetting by row will be effectively ignored, only altering the stored dimensions of \code{x} without changing the values.
Similarly, if the columns are repeated, any subsetting by column will be ignored.
If neither are repeated, subsetting behaves as it would for a normal matrix.

Combining of a CompressedMatrix object will also make use of the repeat structure.
If rows are repeated in all objects to be combined, the output of \code{cbind} will also have repeated rows.
Similarly, if columns are repeated, the output of \code{rbind} will also have repeated columns.
Otherwise, all objects are expanded to their full size prior to combining.

Binary operators work on pairs of CompressedMatrix objects, again preserving the repeat structure whenever possible.
Extracting dimensions uses a second \code{Dims} field in the attributes, bypassing the \code{dim} for a base matrix.
Calling \code{as.matrix} on a \code{CompressedMatrix} object will return the ordinary (uncompressed) matrix.
}

\value{
A object of class CompressedMatrix, containing \code{x} and the additional attributes \code{repeat.row} and \code{repeat.col}.
}

\author{Aaron Lun}

\examples{
# Repeated rows:
library.sizes <- runif(4, 1e6, 2e6)
lib.mat <- makeCompressedMatrix(library.sizes, c(10, 4), byrow=TRUE)
lib.mat

lib.mat[,1:2] # subset by column works as expected
lib.mat[1:10,] # subset by row has no effect (see Details)
as.matrix(lib.mat)

# Repeated columns:
gene.disp <- runif(10, 0.01, 0.1)
disp.mat <- makeCompressedMatrix(gene.disp, c(10, 4), byrow=FALSE)
disp.mat

disp.mat[,1:2] # subset by column has no effect
disp.mat[1:5,] # subset by row works as expected
as.matrix(disp.mat)

# Scalar:
weights <- makeCompressedMatrix(1, c(10, 4))
weights[1:10,] # subsetting has no effect
weights[,1:10]
as.matrix(weights)

# Matrix:
offsets <- makeCompressedMatrix(matrix(runif(40), 10, 4))
offsets[1:5,]
offsets[,1:2]
as.matrix(offsets)
}

\seealso{
\code{\link{as.matrix}}, \code{\link{expandAsMatrix}}
}

\concept{edgeR classes}
