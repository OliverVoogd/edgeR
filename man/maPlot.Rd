\name{maPlot}
\alias{maPlot}
\title{
Plots Log-Fold Change versus Log-Concentration (or, M versus A) for Count Data
}
\description{
To represent counts that were low (e.g. zero in 1 library and non-zero in the other) in one of the two conditions, a 'smear' of points at low A value is presented.
}
\usage{
maPlot(x, y, logAbundance=NULL, logFC=NULL, normalize=FALSE, plot.it=TRUE,
     smearWidth=1, col=NULL, allCol="black", lowCol="orange", deCol="red",
     de.tags=NULL, smooth.scatter=FALSE, lowess=FALSE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of counts or concentrations (group 1)}
  \item{y}{vector of counts or concentrations (group 2)}
  \item{logAbundance}{vector providing the abundance of each gene on the log2 scale. Purely optional (default is \code{NULL}), but in combination with \code{logFC} provides a more direct way to create an MA-plot if the log-abundance and log-fold change are available.}
  \item{logFC}{vector providing the log-fold change for each gene for a given experimental contrast. Default is \code{NULL}, only to be used together with \code{logAbundance} as both need to be non-null for their values to be used.}
  \item{normalize}{logical, whether to divide \code{x} and \code{y} vectors by their sum}
  \item{plot.it}{logical, whether to produce a plot}
  \item{smearWidth}{scalar, width of the smear}
  \item{col}{vector of colours for the points (if \code{NULL}, uses \code{allCol} and \code{lowCol})}
  \item{allCol}{colour of the non-smeared points}
  \item{lowCol}{colour of the smeared points}
  \item{deCol}{colour of the DE (differentially expressed) points}
 \item{de.tags}{indices for genes identified as being differentially expressed; use \code{exactTest} or \code{glmLRT} to identify DE genes. Note that `tag' and `gene' are synonymous here.}
  \item{smooth.scatter}{logical, whether to produce a 'smooth scatter' plot using the graphics::smoothScatter function or just a regular scatter plot; default is \code{FALSE}, i.e. produce a regular scatter plot}
  \item{lowess}{logical, indicating whether or not to add a lowess curve to the MA-plot to give an indication of any trend in the log-fold change with log-concentration}
  \item{\dots}{further arguments passed on to \code{plot}}
}
\details{
The points to be smeared are identified as being equal to the minimum in one of the two groups.  The smear is created by using random uniform numbers of width \code{smearWidth} to the left of the minimum A value.
}
\value{a plot to the current device (if \code{plot.it=TRUE}), and invisibly returns the \code{M} (logFC) and  \code{A} (logConc)  values used for the plot, plus identifiers \code{w} and \code{v} of genes for which \code{M} and {A} values, or just \code{M} values, respectively, were adjusted to make a nicer looking plot. }

\seealso{
\code{\link{plotSmear}}
}

\author{Mark Robinson, Davis McCarthy}

\examples{
y <- matrix(rnbinom(10000,mu=5,size=2),ncol=4)
maPlot(y[,1], y[,2])
}

\concept{Model fit}

