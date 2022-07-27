\name{camera.DGEList}
\alias{camera.DGEList}
\title{Competitive Gene Set Tests for Digital Gene Expression Data}
\description{
Test whether a set of genes is highly ranked relative to other genes in terms of differential expression, accounting for inter-gene correlation.
}

\usage{
\method{camera}{DGEList}(y, index, design, contrast = ncol(design), weights = NULL,
       use.ranks = FALSE, allow.neg.cor=FALSE, inter.gene.cor=0.01, sort = TRUE, \dots)
}

\arguments{
  \item{y}{a \code{DGEList} object containing dispersion estimates.}
  \item{index}{an index vector or a list of index vectors.  Can be any vector such that \code{y[index,]} selects the rows corresponding to the test set.  The list can be made using \code{\link{ids2indices}}.}
  \item{design}{design matrix. Defaults to \code{y$design} or, failing that, to \code{model.matrix(~y$samples$group)}.}
  \item{contrast}{contrast of the linear model coefficients for which the test is required. Can be an integer specifying a column of \code{design}, or else a numeric vector of same length as the number of columns of \code{design}.}
  \item{weights}{numeric matrix of observation weights of same size as \code{y}, or a numeric vector of array weights with length equal to \code{ncol(y)}, or a numeric vector of gene weights with length equal to \code{nrow(y)}.}
  \item{use.ranks}{do a rank-based test (\code{TRUE}) or a parametric test (\code{FALSE})?}
  \item{allow.neg.cor}{should reduced variance inflation factors be allowed for negative correlations?}
  \item{inter.gene.cor}{numeric, optional preset value for the inter-gene correlation within tested sets.  If \code{NA} or \code{NULL}, then an inter-gene correlation will be estimated for each tested set.}
  \item{sort}{logical, should the results be sorted by p-value?}
  \item{\dots}{other arguments are not currently used}
}

\details{
The camera gene set test was proposed by Wu and Smyth (2012) for microarray data.
This function makes the camera test available for digital gene expression data.
The negative binomial count data is converted to approximate normal deviates by computing mid-p quantile residuals (Dunn and Smyth, 1996; Routledge, 1994) under the null hypothesis that the contrast is zero.
See \code{\link{camera}} for more description of the test and for a complete list of possible arguments.
}

\value{
A data.frame giving the gene set results, with most significant sets to the top.
See \code{\link{camera}} for details.
}

\author{Yunshun Chen, Gordon Smyth}

\references{
Dunn, PK, and Smyth, GK (1996).
Randomized quantile residuals.
\emph{J. Comput. Graph. Statist.}, 5, 236-244. 
\url{http://www.statsci.org/smyth/pubs/residual.html}

Routledge, RD (1994).
Practicing safe statistics with the mid-p.
\emph{Canadian Journal of Statistics} 22, 103-110.

Wu, D, and Smyth, GK (2012). Camera: a competitive gene set test accounting for inter-gene correlation.
\emph{Nucleic Acids Research} 40, e133.
\doi{10.1093/nar/gks461}
}

\seealso{
\code{\link{roast.DGEList}},
\code{\link{camera}}.
}

\examples{
mu <- matrix(10, 100, 4)
group <- factor(c(0,0,1,1))
design <- model.matrix(~group)

# First set of 10 genes that are genuinely differentially expressed
iset1 <- 1:10
mu[iset1,3:4] <- mu[iset1,3:4]+10

# Second set of 10 genes are not DE
iset2 <- 11:20

# Generate counts and create a DGEList object
y <- matrix(rnbinom(100*4, mu=mu, size=10),100,4)
y <- DGEList(counts=y, group=group)

# Estimate dispersions
y <- estimateDisp(y, design)

camera(y, iset1, design)
camera(y, iset2, design)

camera(y, list(set1=iset1,set2=iset2), design)
}
