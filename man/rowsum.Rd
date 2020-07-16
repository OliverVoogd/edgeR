\name{rowsum}
\alias{rowsum.DGEList}
\alias{rowsum.SummarizedExperiment}

\title{Sum Over Groups of Genes}
\description{
Condense the rows of a DGEList object so that counts are summed over specified groups of genes.
}
\usage{
\method{rowsum}{DGEList}(x, group, reorder=FALSE, na.rm=FALSE, ...)
\method{rowsum}{SummarizedExperiment}(x, group, reorder=FALSE, na.rm=FALSE, ...)
}
\arguments{
  \item{x}{a \code{DGEList} object or a \code{SummarizedExperiment} object}
  \item{group}{a vector or factor giving the grouping, with one element per row of \code{x}. Missing values will be treated as another group and a warning will be given.}
  \item{reorder}{if \code{TRUE}, then the row.names of the resulting DGEList will be in order of \code{sort(unique(group))}, if \code{FALSE}, they will be in the order that groups were encountered.}
  \item{na.rm}{logical (\code{TRUE} or \code{FALSE}). Should \code{NA} (including \code{NaN}) values be discarded?}
  \item{\dots}{other arguments are not currently used}  
  
}
\details{
If \code{x} is a \code{SummarizedExperiment} object, it is first converted into a \code{DGEList} object. 

A new \code{DGEList} object is computed, with the same columns as \code{x}, but for which the rows correspond to the unique values of \code{group}.
The counts for rows with the same \code{group} value are summed.

Columns of \code{x$genes} will be retained in the output if they contain group-level annotation.
Columns that vary within groups will be dropped.
}
\value{
\code{DGEList} object with the same number of columns as \code{x} and rows corresponding to the unique values of \code{group}.
}
\author{Gordon Smyth}
\seealso{
  \code{\link[base]{rowsum}} in the base package.
}
\examples{
x <- DGEList(counts=matrix(1:8,4,2))
rowsum(x, group=c("A","A","B","B"))
}
