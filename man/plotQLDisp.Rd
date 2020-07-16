\title{Plot the quasi-likelihood dispersion}
\name{plotQLDisp}
\alias{plotQLDisp}
\description{
Plot the genewise quasi-likelihood dispersion against the gene abundance (in log2 counts per million).
}
\usage{
plotQLDisp(glmfit, xlab="Average Log2 CPM", ylab="Quarter-Root Mean Deviance", pch=16, 
       cex=0.2, col.shrunk="red", col.trend="blue", col.raw="black", \dots)
}
\arguments{
  \item{glmfit}{a \code{DGEGLM} object produced by \code{\link{glmQLFit}}.}
  \item{xlab}{label for the x-axis.}
  \item{ylab}{label for the y-axis.}
  \item{pch}{the plotting symbol. See \code{\link{points}} for more details.}
  \item{cex}{plot symbol expansion factor. See \code{\link{points}} for more details.}
  \item{col.shrunk}{color of the points representing the squeezed quasi-likelihood dispersions.}
  \item{col.trend}{color of line showing dispersion trend.}
  \item{col.raw}{color of points showing the unshrunk dispersions.}
  \item{\dots}{any other arguments are passed to \code{plot}.}
}

\details{
This function displays the quarter-root of the quasi-likelihood dispersions for all genes, before and after shrinkage towards a trend.
If \code{glmfit} was constructed without an abundance trend, the function instead plots a horizontal line (of colour \code{col.trend}) at the common value towards which dispersions are shrunk.
The quarter-root transformation is applied to improve visibility for dispersions around unity.
}

\value{
A plot is created on the current graphics device.
}

\author{Aaron Lun, Davis McCarthy, Gordon Smyth, Yunshun Chen.}

\examples{
nbdisp <- 1/rchisq(1000, df=10)
y <- DGEList(matrix(rnbinom(6000, size = 1/nbdisp, mu = 10),1000,6))
design <- model.matrix(~factor(c(1,1,1,2,2,2)))
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)
plotQLDisp(fit)

fit <- glmQLFit(y, design, abundance.trend=FALSE)
plotQLDisp(fit)
}

\references{
Chen Y, Lun ATL, and Smyth, GK (2016).
From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline.
\emph{F1000Research} 5, 1438.
\url{http://f1000research.com/articles/5-1438}
}

\concept{Dispersion estimation}
