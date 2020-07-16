\name{subsetting}
\alias{subsetting}
\alias{[.DGEList}
\alias{[.DGEExact}
\alias{[.DGELRT}
\alias{[.DGEGLM}
\alias{[.TopTags}
\title{Subset DGEList, DGEGLM, DGEExact and DGELRT Objects}
\description{
Extract a subset of a \code{DGEList}, \code{DGEGLM}, \code{DGEExact} or \code{DGELRT} object.
}
\usage{
\method{[}{DGEList}(object, i, j, keep.lib.sizes=TRUE)
\method{[}{DGEGLM}(object, i, j)
\method{[}{DGEExact}(object, i, j)
\method{[}{DGELRT}(object, i, j)
\method{[}{TopTags}(object, i, j)
}
\arguments{
  \item{object}{object of class \code{DGEList}, \code{DGEGLM}, \code{DGEExact} or \code{DGELRT}. For \code{subsetListOfArrays}, any list of conformal matrices and vectors.}
  \item{i,j}{elements to extract. \code{i} subsets the genes while \code{j} subsets the libraries.
  Note that columns of \code{DGEGLM}, \code{DGEExact} and \code{DGELRT} objects cannot be subsetted.}
  \item{keep.lib.sizes}{logical, if \code{TRUE} the lib.sizes will be kept unchanged on output, otherwise they will be recomputed as the column sums of the counts of the remaining rows.}
}
\details{
\code{i,j} may take any values acceptable for the matrix components of \code{object} of class \code{DGEList}.
See the \link{Extract} help entry for more details on subsetting matrices. For \code{DGEGLM}, \code{DGEExact} and \code{DGELRT} objects, only rows (i.e. \code{i}) may be subsetted.
}
\value{
An object of the same class as \code{object} holding data from the specified subset of rows and columns.
}
\author{Davis McCarthy, Gordon Smyth}
\seealso{
  \code{\link{Extract}} in the base package.
}
\examples{
d <- matrix(rnbinom(16,size=1,mu=10),4,4)
rownames(d) <- c("a","b","c","d")
colnames(d) <- c("A1","A2","B1","B2")
d <- DGEList(counts=d,group=factor(c("A","A","B","B")))
d[1:2,]
d[1:2,2]
d[,2]
d <- estimateCommonDisp(d)
results <- exactTest(d)
results[1:2,]
# NB: cannot subset columns for DGEExact objects
}

\concept{edgeR classes}

