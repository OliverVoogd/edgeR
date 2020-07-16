\name{splitIntoGroups}
\alias{splitIntoGroups}
\alias{splitIntoGroups.DGEList}
\alias{splitIntoGroups.default}
\alias{splitIntoGroupsPseudo}

\title{Split the Counts or Pseudocounts from a DGEList Object According To Group}

\description{Split the counts from a DGEList object according to group, creating a list where each element consists of a numeric matrix of counts for a particular experimental group. Given a pair of groups, split pseudocounts for these groups, creating a list where each element is a matrix of pseudocounts for a particular gourp.}

\usage{
\S3method{splitIntoGroups}{DGEList}(y, ...)
\S3method{splitIntoGroups}{default}(y, group=NULL, ...)
splitIntoGroupsPseudo(pseudo, group, pair)
}

\arguments{
\item{y}{matrix of counts or a \code{DGEList} object.}
\item{group}{vector or factor giving the experimental group/condition for each library.}
\item{pseudo}{numeric matrix of quantile-adjusted pseudocounts to be split}
\item{pair}{vector of length two stating pair of groups to be split for the pseudocounts}
\item{\dots}{other arguments that are not currently used.}
}

\value{\code{splitIntoGroups} outputs a list in which each element is a matrix of count counts for an individual group. \code{splitIntoGroupsPseudo} outputs a list with two elements, in which each element is a numeric matrix of (pseudo-)count data for one of the groups specified.}

\author{Davis McCarthy}
\examples{
# generate raw counts from NB, create list object
y <- matrix(rnbinom(80, size=1, mu=10), nrow=20)
d <- DGEList(counts=y, group=rep(1:2, each=2), lib.size=rep(c(1000:1001), 2))
rownames(d$counts) <- paste("gene", 1:nrow(d$counts), sep=".")
z1 <- splitIntoGroups(d)
z2 <- splitIntoGroupsPseudo(d$counts, d$group, pair=c(1,2))
}

