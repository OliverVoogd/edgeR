\name{diffSpliceDGE}
\alias{diffSpliceDGE}

\title{Test for Differential Exon Usage}
\description{Given a negative binomial generalized log-linear model fit at the exon level, test for differential exon usage between experimental conditions.}
\usage{
diffSpliceDGE(glmfit, coef=ncol(glmfit$design), contrast=NULL, geneid, exonid=NULL,
              prior.count=0.125, verbose=TRUE)
}

\arguments{
  \item{glmfit}{an \code{DGEGLM} fitted model object produced by \code{glmFit} or \code{glmQLFit}. Rows should correspond to exons.}
  \item{coef}{integer indicating which coefficient of the generalized linear model is to be tested for differential exon usage. Defaults to the last coefficient.}
  \item{contrast}{numeric vector specifying the contrast of the linear model coefficients to be tested for differential exon usage. Length must equal to the number of columns of \code{design}. If specified, then takes precedence over \code{coef}.}
  \item{geneid}{gene identifiers. Either a vector of length \code{nrow(glmfit)} or the name of the column of \code{glmfit$genes} containing the gene identifiers. Rows with the same ID are assumed to belong to the same gene.}
  \item{exonid}{exon identifiers. Either a vector of length \code{nrow(glmfit)} or the name of the column of \code{glmfit$genes} containing the exon identifiers.}
  \item{prior.count}{average prior count to be added to observation to shrink the estimated log-fold-changes towards zero.}
  \item{verbose}{logical, if \code{TRUE} some diagnostic information about the number of genes and exons is output.}
}

\value{
\code{diffSpliceDGE} produces an object of class \code{DGELRT} containing the component \code{design} from \code{glmfit} plus the following new components:
  \item{comparison}{character string describing the coefficient being tested.}
  \item{coefficients}{numeric vector of coefficients on the natural log scale. Each coefficient is the difference between the log-fold-change for that exon versus the log-fold-change for the entire gene which contains that exon.}
  \item{genes}{data.frame of exon annotation.}
  \item{genecolname}{character string giving the name of the column of \code{genes} containing gene IDs.}
  \item{exoncolname}{character string giving the name of the column of \code{genes} containing exon IDs.}
  \item{exon.df.test}{numeric vector of testing degrees of freedom for exons.}
  \item{exon.p.value}{numeric vector of p-values for exons.}
  \item{gene.df.test}{numeric vector of testing degrees of freedom for genes.}
  \item{gene.p.value}{numeric vector of gene-level testing p-values.}
  \item{gene.Simes.p.value}{numeric vector of Simes' p-values for genes.}
  \item{gene.genes}{data.frame of gene annotation.}

Some components of the output depend on whether \code{glmfit} is produced by \code{glmFit} or \code{glmQLFit}. 
If \code{glmfit} is produced by \code{glmFit}, then the following components are returned in the output object:
  \item{exon.LR}{numeric vector of LR-statistics for exons.}
  \item{gene.LR}{numeric vector of LR-statistics for gene-level test.}
  
If \code{glmfit} is produced by \code{glmQLFit}, then the following components are returned in the output object:
  \item{exon.F}{numeric vector of F-statistics for exons.}
  \item{gene.df.prior}{numeric vector of prior degrees of freedom for genes.}
  \item{gene.df.residual}{numeric vector of residual degrees of freedom for genes.}
  \item{gene.F}{numeric vector of F-statistics for gene-level test.}

The information and testing results for both exons and genes are sorted by geneid and by exonid within gene.
}

\details{
This function tests for differential exon usage for each gene for a given coefficient of the generalized linear model.

Testing for differential exon usage is equivalent to testing whether the exons in each gene have the same log-fold-changes as the other exons in the same gene. 
At exon-level, the log-fold-change of each exon is compared to the log-fold-change of the entire gene which contains that exon.
At gene-level, two different tests are provided. One is converting exon-level p-values to gene-level p-values by the Simes method.
The other is using exon-level test statistics to conduct gene-level tests.
}

\author{Yunshun Chen and Gordon Smyth}

\examples{
# Gene exon annotation
Gene <- paste("Gene", 1:100, sep="")
Gene <- rep(Gene, each=10)
Exon <- paste("Ex", 1:10, sep="")
Gene.Exon <- paste(Gene, Exon, sep=".")
genes <- data.frame(GeneID=Gene, Gene.Exon=Gene.Exon)

group <- factor(rep(1:2, each=3))
design <- model.matrix(~group)
mu <- matrix(100, nrow=1000, ncol=6)
# knock-out the first exon of Gene1 by 90%
mu[1,4:6] <- 10
# generate exon counts
counts <- matrix(rnbinom(6000,mu=mu,size=20),1000,6)

y <- DGEList(counts=counts, lib.size=rep(1e6,6), genes=genes)
gfit <- glmFit(y, design, dispersion=0.05)

ds <- diffSpliceDGE(gfit, geneid="GeneID")
topSpliceDGE(ds)
plotSpliceDGE(ds)
}

\concept{Differential exon usage}
