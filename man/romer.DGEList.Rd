\name{romer.DGEList}
\alias{romer.DGEList}
\title{Rotation Gene Set Enrichment for Digital Gene Expression Data}
\description{
Romer gene set enrichment tests for Negative Binomial generalized linear models.
}

\usage{
\method{romer}{DGEList}(y, index, design=NULL, contrast=ncol(design), \dots)
}

\arguments{
  \item{y}{\code{DGEList} object.}
  \item{index}{list of indices specifying the rows of \code{y} in the gene sets. The list can be made using \link{ids2indices}.}
  \item{design}{design matrix. Defaults to \code{y$design} or, failing that, to \code{model.matrix(~y$samples$group)}.}
  \item{contrast}{contrast for which the test is required. Can be an integer specifying a column of \code{design}, or the name of a column of \code{design}, or else a contrast vector of length equal to the number of columns of \code{design}.}
  \item{\dots}{other arguments are passed to \code{\link{romer.default}}. For example, the number of rotations \code{nrot} can be increased from the default of \code{9999} to increase the resolution of the p-values.}
}

\value{
Numeric matrix giving p-values and the number of matched genes in each gene set.
Rows correspond to gene sets.
There are four columns giving the number of genes in the set and p-values for the alternative hypotheses up, down or mixed.
See \code{\link{romer}} for details.
}

\details{
The ROMER procedure described by Majewski et al (2010) is implemented in \code{romer} in the limma package.
This \code{romer} method for \code{DGEList} objects makes the romer procedure available for count data such as RNA-seq data.
The negative binomial count data is converted to approximate normal deviates by computing mid-p quantile residuals (Dunn and Smyth, 1996; Routledge, 1994) under the null hypothesis that the contrast is zero.
The normal deviates are then passed to the \code{romer} function in limma.
See \code{\link{romer}} for more description of the test and for a complete list of possible arguments.
}

\seealso{
\code{\link{romer}}
}

\author{Yunshun Chen and Gordon Smyth}

\references{
Majewski, IJ, Ritchie, ME, Phipson, B, Corbin, J, Pakusch, M, Ebert, A, Busslinger, M, Koseki, H, Hu, Y, Smyth, GK, Alexander, WS, Hilton, DJ, and Blewitt, ME (2010).
Opposing roles of polycomb repressive complexes in hematopoietic stem and progenitor cells.
\emph{Blood}, 116, 731-719.
\url{http://www.ncbi.nlm.nih.gov/pubmed/20445021}

Dunn, PK, and Smyth, GK (1996).
Randomized quantile residuals.
\emph{J. Comput. Graph. Statist.}, 5, 236-244. 
\url{http://www.statsci.org/smyth/pubs/residual.html}

Routledge, RD (1994).
Practicing safe statistics with the mid-p.
\emph{Canadian Journal of Statistics} 22, 103-110.
}

\examples{
mu <- matrix(10, 100, 4)
group <- factor(c(0,0,1,1))
design <- model.matrix(~group)

# First set of 10 genes that are genuinely differentially expressed
iset1 <- 1:10
mu[iset1,3:4] <- mu[iset1,3:4]+20

# Second set of 10 genes are not DE
iset2 <- 11:20

# Generate counts and create a DGEList object
y <- matrix(rnbinom(100*4, mu=mu, size=10),100,4)
y <- DGEList(counts=y, group=group)

# Estimate dispersions
y <- estimateDisp(y, design)

romer(y, iset1, design, contrast=2)
romer(y, iset2, design, contrast=2)
romer(y, list(set1=iset1, set2=iset2), design, contrast=2)
}

\concept{Gene set testing}
