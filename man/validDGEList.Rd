\name{validDGEList}
\alias{validDGEList}
\alias{validDGEList}

\title{Check for Valid DGEList object}

\description{
Check for existence of standard components of DGEList object.
}

\usage{
validDGEList(y)
}

\arguments{
\item{y}{\code{DGEList} object.}
}

\details{
This function checks that the standard \code{counts} and \code{samples} components of a \code{DGEList} object are present.
}

\value{
\code{DGEList} with missing components added. 
}

\author{Gordon Smyth}

\examples{
counts <- matrix(rpois(4*2,lambda=5),4,2)
dge <- new("DGEList", list(counts=counts))
validDGEList(dge)
}

\seealso{
\code{\link{DGEList}}
}

\concept{edgeR classes}
