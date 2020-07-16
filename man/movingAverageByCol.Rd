\name{movingAverageByCol}
\alias{movingAverageByCol}
\title{Moving Average Smoother of Matrix Columns}
\description{
Apply a moving average smoother to the columns of a matrix.
}
\usage{
movingAverageByCol(x, width=5, full.length=TRUE)
}
\arguments{
  \item{x}{numeric matrix}
  \item{width}{integer, width of window of rows to be averaged}
  \item{full.length}{logical value, should output have same number of rows as input?}
}
\details{
If \code{full.length=TRUE}, narrower windows are used at the start and end of each column to make a column of the same length as input.
If \code{FALSE}, all values are averager of \code{width} input values, so the number of rows is less than input.
}
\value{
Numeric matrix containing smoothed values.
If \code{full.length=TRUE}, of same dimension as \code{x}.
If \code{full.length=FALSE}, has \code{width-1} fewer rows than \code{x}.
}
\examples{
x <- matrix(rpois(20,lambda=5),10,2)
movingAverageByCol(x,3)
}
\author{Gordon Smyth}

