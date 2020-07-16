\name{head}
\alias{head.DGEList}
\alias{head.DGEExact}
\alias{head.DGEGLM}
\alias{head.DGELRT}
\alias{head.TopTags}
\alias{tail.DGEList}
\alias{tail.DGEExact}
\alias{tail.DGEGLM}
\alias{tail.DGELRT}
\alias{tail.TopTags}

\title{Return the First to Last Part of a Data Object}

\description{
Retrieve the first or last parts of a DGEList, DGEExat, DGEGLM, DGELRT or TopTags object.
}

\usage{
\method{head}{DGEList}(x, n = 6L, \dots)
\method{tail}{DGEList}(x, n = 6L, \dots)
}

\arguments{
  \item{x}{an object of class \code{DGEList}, \code{DGEExact}, \code{DGEGLM}, \code{DGELRT} or \code{TopTags}.}
  \item{n}{
    a single integer.
    If positive or zero, number rows of resulting object.
    If negative, all but the \code{n} last/first rows of \code{x}.
  }
  \item{\dots}{other arguments are not currently used.}
}

\details{
\code{head} (\code{tail}) returns the first (last) \code{n} rows when \code{n >= 0} or all but the last (first) \code{n} rows when \code{n < 0}.
}

\value{
An object like \code{x} but generally with fewer rows.
}

\author{Gordon Smyth}

\seealso{
  \code{\link{head}} in the utils package or \code{\link{head.EList}} in the limma package.
}

\examples{
Counts <- matrix(rpois(40,lambda=10),20,2)
rownames(Counts) <- paste0("Gene",1:20)
colnames(Counts) <- c("A","B")
y <- DGEList(Counts)
head(y)
tail(y)
}

\concept{edgeR classes}
