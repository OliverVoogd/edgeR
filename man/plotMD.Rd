\title{Mean-Difference Plot of Count Data}
\name{plotMD.DGEList}
\alias{plotMD.DGEList}
\alias{plotMD.SummarizedExperiment}
\alias{plotMD.DGEGLM}
\alias{plotMD.DGELRT}
\alias{plotMD.DGEExact}

\description{
Creates a mean-difference plot (aka MA plot) with color coding for highlighted points.
}

\usage{
\method{plotMD}{DGEList}(object, column=1, xlab="Average log CPM (this sample and others)",
       ylab="log-ratio (this sample vs others)",
       main=colnames(object)[column], status=object$genes$Status,
       zero.weights=FALSE, prior.count=3, \dots)
\method{plotMD}{SummarizedExperiment}(object, column=1, xlab="Average log CPM (this sample and others)", 
       ylab="log-ratio (this sample vs others)", zero.weights=FALSE, prior.count=3, ...)
\method{plotMD}{DGEGLM}(object, column=ncol(object), coef=NULL, xlab="Average log CPM",
       ylab="log-fold-change", main=colnames(object)[column],
       status=object$genes$Status, zero.weights=FALSE, \dots)
\method{plotMD}{DGELRT}(object, xlab="Average log CPM", ylab="log-fold-change", 
       main=object$comparison, status=object$genes$Status, contrast=1, 
       adjust.method="BH", p.value=0.05, \dots)
\method{plotMD}{DGEExact}(object, xlab="Average log CPM", ylab="log-fold-change", 
       main=NULL, status=object$genes$Status,
       adjust.method="BH", p.value=0.05, \dots)
}

\arguments{
  \item{object}{an object of class \code{DGEList}, \code{SummarizedExperiment}, \code{DGEGLM}, \code{DGEGLM} or \code{DGEExact}.}
  \item{column}{integer, column of \code{object} to be plotted.}
  \item{coef}{alternative to \code{column} for fitted model objects. If specified, then \code{column} is ignored.}
  \item{xlab}{character string, label for x-axis}
  \item{ylab}{character string, label for y-axis}
  \item{main}{character string, title for plot}
  \item{status}{vector giving the control status of each spot on the array, of same length as the number of rows of \code{object}.
  If \code{NULL} under the \code{DGEList}, \code{SummarizedExperiment} or \code{DGEGLM} method, then all points are plotted in the default color, symbol and size.
  If \code{NULL} under the \code{DGELRT} or \code{DGEExact} method, then \code{\link{decideTestsDGE}} is run to determine the status of all the genes. The up-regulated  DE genes are highlighted in red and down-regulated in blue.}
  \item{zero.weights}{logical, should spots with zero or negative weights be plotted?}
  \item{prior.count}{the average prior count to be added to each observation. Larger values produce more shrinkage.}
  \item{contrast}{integer specifying which log-fold-change to be plotted in the case of testing multiple contrasts. Only used for the \code{DGELRT} method with multiple contrasts.}
  \item{adjust.method}{character string passed to \code{\link{decideTestsDGE}} specifying p-value adjustment method. Only used when \code{status} is \code{NULL}. See \code{\link{decideTestsDGE}} for details.}
  \item{p.value}{numeric value between 0 and 1 giving the desired size of the test. Only used and passed to \code{\link{decideTestsDGE}} when \code{status} is \code{NULL}.}
  \item{\dots}{other arguments are passed to \code{\link{plotWithHighlights}}.}
}

\details{
A mean-difference plot (MD-plot) is a plot of log fold changes (differences) versus average log values (means).
The history of mean-difference plots and MA-plots is reviewed in Ritchie et al (2015).

For \code{DGEList} and \code{SummarizedExperiment} objects, a between-sample MD-plot is produced.
Counts are first converted to log2-CPM values.
An articifial array is produced by averaging all the samples other than the sample specified.
A mean-difference plot is then producing from the specified sample and the artificial sample.
This procedure reduces to an ordinary mean-difference plot when there are just two arrays total.

If \code{object} is an \code{DGEGLM} object, then the plot is an fitted model MD-plot in which the estimated coefficient is on the y-axis and the average logCPM value is on the x-axis.
If \code{object} is an \code{DGEExact} or \code{DGELRT} object, then the MD-plot displays the logFC vs the logCPM values from the results table.

The \code{status} vector can correspond to any grouping of the probes that is of interest.
If \code{object} is a fitted model object, then \code{status} vector is often used to indicate statistically significance, so that differentially expressed points are highlighted.

The \code{status} can be included as the component \code{object$genes$Status} instead of being passed as an argument to \code{plotMD}.

See \code{\link{plotWithHighlights}} for how to set colors and graphics parameters for the highlighted and non-highlighted points.
}

\value{A plot is created on the current graphics device.}

\references{
Ritchie, ME, Phipson, B, Wu, D, Hu, Y, Law, CW, Shi, W, and Smyth, GK (2015).
limma powers differential expression analyses for RNA-sequencing and microarray studies.
\emph{Nucleic Acids Research} Volume 43, e47.
\doi{10.1093/nar/gkv007}
}

\author{Gordon Smyth}

\seealso{
\code{plotSmear}

The driver function for \code{plotMD} is \code{\link{plotWithHighlights}}.
}

\concept{Model fit}
