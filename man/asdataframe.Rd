\name{as.data.frame}
\alias{as.data.frame.DGEExact}
\alias{as.data.frame.DGELRT}
\alias{as.data.frame.TopTags}
\title{Turn a TopTags Object into a Dataframe}
\description{
Turn a \code{TopTags} object into a \code{data.frame}.
}
\usage{
\method{as.data.frame}{TopTags}(x, row.names = NULL, \dots)
}
\arguments{
  \item{x}{an object of class \code{TopTags}}
  \item{row.names}{\code{NULL} or a character vector giving the row names for the
          data frame.  Missing values are not allowed.}
  \item{\dots}{other arguments are not currently used.}
}
\details{
Convert edgeR objects into data.frames.
This method returns the \code{table} component of a \code{TopTags} object.
For \code{DGEExact} and \code{DGELRT} objects, the \code{genes} and \code{table} components are combined into a data.frame, similar to what is done by \code{topTags} but without sorting or p-value adjustment.
}
\value{
A data.frame.
}
\author{Gordon Smyth}
\seealso{
  \code{\link[base]{as.data.frame}} in the base package.
}

\concept{edgeR classes}
