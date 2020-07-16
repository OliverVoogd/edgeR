\name{cbind}
\alias{cbind}
\alias{cbind.DGEList}
\alias{rbind.DGEList}
\title{Combine DGEList Objects}
\description{
Combine a set of \code{DGEList} objects.
}
\usage{
\method{cbind}{DGEList}(\dots, deparse.level=1)
\method{rbind}{DGEList}(\dots, deparse.level=1)
}
\arguments{
  \item{\dots}{\code{DGEList} objects.}
  \item{deparse.level}{not currently used, see \code{\link[base]{cbind}} in the base package}
}
\details{
\code{cbind} combines data objects assuming the same genes in the same order but different samples.
\code{rbind} combines data objects assuming equivalent samples, i.e., the same RNA targets, but different genes.

For \code{cbind}, the matrices of count data from the individual objects are cbinded.
The data.frames of samples information, if they exist, are rbinded.
The combined data object will preserve any additional components or attributes found in the first object to be combined.
For \code{rbind}, the matrices of count data are rbinded while the sample information is unchanged.
}
\value{
An \code{\link[edgeR:DGEList]{DGEList}}
object holding data from all samples and all genes from the individual objects.
}
\author{Gordon Smyth}
\seealso{
  \code{\link[base]{cbind}} in the base package.
}
\examples{
\dontrun{
dge <- cbind(dge1,dge2,dge3)
}
}

\concept{edgeR classes}
