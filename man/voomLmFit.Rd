\name{voomLmFit}
\alias{voomLmFit}
\title{Apply voom-lmFit Pipeline While Accounting for Loss of Residual DF Due to Exact Zeros}
\description{
Transform count data to log2-counts per million (logCPM), estimate voom precision weights and fit limma linear models while allowing for loss of residual degrees of freedom due to exact zeros.
}

\usage{
voomLmFit(counts, design = NULL, block = NULL, prior.weights = NULL,
     sample.weights = FALSE, var.design = NULL, var.group = NULL, 
     lib.size = NULL, normalize.method = "none",
     span = 0.5, plot = FALSE, save.plot = FALSE, keep.EList = TRUE)
}
\arguments{
  \item{counts}{
    a numeric \code{matrix} containing raw counts, or an \code{ExpressionSet} containing raw counts, or a \code{DGEList} object.
    Counts must be non-negative and NAs are not permitted.}
  \item{design}{
    design matrix with rows corresponding to samples and columns to coefficients to be estimated.
    Defaults to the unit vector meaning that samples are treated as replicates.}
  \item{block}{
    vector or factor specifying a blocking variable on the samples.
    Has length equal to the number of \code{ncol(counts)}.
    Samples within each block are assumed to be correlated.}
  \item{prior.weights}{
    prior weights.
    Can be a numeric matrix of individual weights of same dimensions as the \code{counts},
    or a numeric vector of sample weights with length equal to \code{ncol(counts)},
    or a numeric vector of gene weights with length equal to \code{nrow(counts)}.}
  \item{sample.weights}{
    logical value, if \code{TRUE} then empirical sample quality weights will be estimated.}
  \item{var.design}{
    optional design matrix for the sample weights.
    Defaults to the sample-specific model whereby each sample has a distinct variance.}
  \item{var.group}{
    optional vector or factor indicating groups to have different array weights.
    This is another way to specify \code{var.design} for groupwise sample weights.}
  \item{lib.size}{
    numeric vector containing total library sizes for each sample.
    Defaults to the normalized (effective) library sizes in \code{counts} if \code{counts} is a \code{DGEList} or to the columnwise count totals if \code{counts} is a matrix.}
  \item{normalize.method}{
    the microarray-style normalization method to be applied to the logCPM values (if any).
    Choices are as for the \code{method} argument of \code{normalizeBetweenArrays} when the data is single-channel.
    Any normalization factors found in \code{counts} will still be used even if \code{normalize.method="none"}.}
  \item{span}{
    width of the smoothing window used for the lowess mean-variance trend.
    Expressed as a proportion between 0 and 1.}
  \item{plot}{
    logical, should a plot of the mean-variance trend be displayed?}
  \item{save.plot}{
    logical, should the coordinates and line of the plot be saved in the output?}
  \item{keep.EList}{
    logical. If \code{TRUE}, then the normalized log2-CPM values and voom weights will be saved in the component \code{EList} of the output object.}
}

\details{
This function adapts the limma voom method (Law et al, 2014) to allow for loss of residual degrees of freedom due to exact zero counts (Lun and Smyth, 2017).
The loss residual df occurs when all the counts in a group are zero or when there are blocking factors that can fit zero counts exactly.
The function transforms the counts to the log2-CPM scale, computes voom precision weights and fits limma linear models.
Residual df are computed similarly as far \code{\link{glmQLFit}}.

The function is analogous to calling \code{voom} followed by \code{duplicateCorrelation} and \code{lmFit} except for the modified residual df values and residual standard deviation \code{sigma} values.
This function returns \code{df.residual} values that are less than or equal to those from \code{lmFit} and \code{sigma} values that are greater than or equal to those from \code{lmFit}.
\code{voomLmFit} is more robust to zero counts than calling \code{voom}, \code{duplicateCorrelation} and \code{lmFit} separately and provides more rigorous error rate control.

If \code{block} is specified, then the intra-block correlation is estimated using \code{duplicateCorrelation}
In that case, the voom weights and the intra-block correlation are each estimated twice to achieve effective convergence.

Empirical sample quality weights will be estimated if \code{sample.weights=TRUE} or if \code{var.design} or \code{var.group} are non-NULL (Liu et al 2015).
In that case, \code{voomLmFit} is analogous to running \code{voomWithQualityWeights} followed by \code{lmFit}.

\code{voomLmFit} is usually followed by running \code{\link{eBayes}} on the fitted model object.
}

\value{
An MArrayLM object containing linear model fits for each row of data.
The object includes a \code{targets} data.frame component containing sample annotation.
Columns of \code{targets} include \code{lib.size} and \code{sample.weight} (if \code{sample.weights=TRUE}).

If \code{save.plot=TRUE} then the output object will include components \code{voom.xy} and \code{voom.line}.
\code{voom.xy} contains the x and y coordinates of the points in the voom mean-variance plot in the same format as produced by \code{\link{xy.coords}} and \code{voom.line} contained the estimated trend curve.

If \code{keep.EList=TRUE} then the output object includes component \code{EList}, which is an \code{EList} object in the same format as produced by \code{voom} containing the voom log2-CPM values and the voom weights.
}

\author{Gordon Smyth}

\references{
Law, CW, Chen, Y, Shi, W, Smyth, GK (2014).
Voom: precision weights unlock linear model analysis tools for RNA-seq read counts.
\emph{Genome Biology} 15, R29.
See also the Preprint Version at \url{http://www.statsci.org/smyth/pubs/VoomPreprint.pdf} incorporating some notational corrections.

Lun, ATL, Smyth, GK (2017).
No counts, no variance: allowing for loss of degrees of freedom when assessing biological variability from RNA-seq data.
\emph{Statistical Applications in Genetics and Molecular Biology} 16(2), 83-93.
\doi{10.1515/sagmb-2017-0010}

Liu, R, Holik, AZ, Su, S, Jansz, N, Chen, K, Leong, HS, Blewitt, ME, Asselin-Labat, M-L, Smyth, GK, Ritchie, ME (2015).
Why weight? Modelling sample and observational level variability improves power in RNA-seq analyses.
\emph{Nucleic Acids Research} 43, e97.
\doi{10.1093/nar/gkv412}
}

\seealso{
\code{\link{voom}},
\code{\link{lmFit}},
\code{\link{voomWithQualityWeights}},
\code{\link{duplicateCorrelation}},
\code{\link{arrayWeights}},
\code{\link{MArrayLM-class}}.
}

\concept{Model fit}
