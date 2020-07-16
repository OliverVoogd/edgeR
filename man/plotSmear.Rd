\name{plotSmear}
\alias{plotSmear}
\title{
Smear plot
}
\description{
Make a mean-difference plot of two libraries of count data with smearing of points with very low counts, especially those that are zero for one of the columns.
}

\usage{
plotSmear(object, pair=NULL, de.tags=NULL, xlab="Average logCPM", ylab="logFC", pch=19,
     cex=0.2, smearWidth=0.5, panel.first=grid(), smooth.scatter=FALSE, lowess=FALSE, \dots)
}
\arguments{
 \item{object}{\code{DGEList}, \code{DGEExact} or \code{DGELRT} object containing data to produce an MA-plot.}
  \item{pair}{pair of experimental conditions to plot (if \code{NULL}, the first two conditions are used). Ignored if \code{object} is a \code{DGELRT} object.}
  \item{de.tags}{rownames for genes identified as being differentially expressed; use \code{exactTest} or \code{glmLRT} to identify DE genes. Note that `tag' and `gene' are synonymous here.}
  \item{xlab}{x-label of plot}
  \item{ylab}{y-label of plot}
  \item{pch}{scalar or vector giving the character(s) to be used in the plot; default value of \code{19} gives a round point.}
  \item{cex}{character expansion factor, numerical value giving the amount by which plotting text and symbols should be magnified relative to the default; default \code{cex=0.2} to make the plotted points smaller}
  \item{smearWidth}{width of the smear}
  \item{panel.first}{an expression to be evaluated after the plot axes are set up but before any plotting takes place; the default \code{grid()} draws a background grid to aid interpretation of the plot}
  \item{smooth.scatter}{logical, whether to produce a 'smooth scatter' plot using the \code{KernSmooth::smoothScatter} function or just a regular scatter plot; default is \code{FALSE}, i.e. produce a regular scatter plot}
  \item{lowess}{logical, indicating whether or not to add a lowess curve to the MA-plot to give an indication of any trend in the log-fold change with log-concentration}
  \item{\dots}{further arguments passed on to \code{plot}}
}

\value{
Invisibly returns the x and y coordinates of the plotted points, and a plot is created on the current device.
}

\details{
\code{plotSmear} produces a type of mean-difference plot (or MA plot) with a special representation (smearing) of log-ratios that are infinite.
\code{plotSmear} resolves the problem of plotting genes that have a total count of zero for one of the groups by adding the 'smear' of points at low A value.
The points to be smeared are identified as being equal to the minimum estimated concentration in one of the two groups.
The smear is created by using random uniform numbers of width \code{smearWidth} to the left of the minimum A.
\code{plotSmear} also allows easy highlighting of differentially expressed (DE) genes.
}

\author{Mark Robinson created the original concept of smearing the infinite log-fold-changes.}

\seealso{
\code{\link{maPlot}}, \code{\link{plotMD.DGEList}}
}

\examples{
y <- matrix(rnbinom(10000,mu=5,size=2),ncol=4)
d <- DGEList(counts=y, group=rep(1:2,each=2), lib.size=colSums(y))
rownames(d$counts) <- paste("gene",1:nrow(d$counts),sep=".")
d <- estimateCommonDisp(d)
plotSmear(d)

# find differential expression
de <- exactTest(d)

# highlighting the top 500 most DE genes
de.genes <- rownames(topTags(de, n=500)$table)
plotSmear(d, de.tags=de.genes)
}

\concept{Data exploration}
\concept{Model fit}
