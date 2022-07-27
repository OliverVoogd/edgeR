\name{DGEList}
\alias{DGEList}

\title{
DGEList Constructor
}

\description{
Creates a \code{DGEList} object from a table of counts (rows=features, columns=samples), group indicator for each column, library size (optional) and a table of feature annotation (optional).
}

\usage{
DGEList(counts = matrix(0, 0, 0), lib.size = colSums(counts),
        norm.factors = rep(1,ncol(counts)), samples = NULL,
        group = NULL, genes = NULL, remove.zeros = FALSE)
}

\arguments{
  \item{counts}{numeric matrix of read counts.}
  \item{lib.size}{numeric vector giving the total count (sequence depth) for each library.}
  \item{norm.factors}{numeric vector of normalization factors that modify the library sizes.}
  \item{samples}{data frame containing information for each sample.}
  \item{group}{vector or factor giving the experimental group/condition for each sample/library.}
  \item{genes}{data frame containing annotation information for each gene.}
  \item{remove.zeros}{logical, whether to remove rows that have 0 total count.}
}

\details{
To facilitate programming pipelines,
\code{NULL} values can be input for \code{lib.size}, \code{norm.factors}, \code{samples} or \code{group}, in which case the default value is used as if the argument had been missing.
}

\value{a \code{\link[edgeR:DGEList-class]{DGEList}} object}

\author{edgeR team. First created by Mark Robinson.}

\seealso{\code{\link[edgeR:DGEList-class]{DGEList-class}}}

\examples{
ngenes <- 1000
nsamples <- 4
Counts <- matrix(rnbinom(ngenes*nsamples,mu=5,size=2),ngenes,nsamples)
rownames(Counts) <- 1:ngenes
y <- DGEList(counts=Counts, group=rep(1:2,each=2))
dim(y)
colnames(y)
y$samples
y$genes <- data.frame(Symbol=paste0("Gene",1:ngenes))
show(y)
}

\concept{edgeR classes}
