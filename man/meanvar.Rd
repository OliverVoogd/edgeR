\name{meanvar}
\alias{binMeanVar}
\alias{plotMeanVar}

\title{Explore the Mean-Variance Relationship for DGE Data}

\description{Appropriate modelling of the mean-variance relationship in DGE data is important for making inferences about differential expression. Here are functions to compute gene means and variances, as well at looking at these quantities when data is binned based on overall expression level.}


\usage{
plotMeanVar(object, meanvar=NULL, show.raw.vars=FALSE, show.tagwise.vars=FALSE,
            show.binned.common.disp.vars=FALSE, show.ave.raw.vars=TRUE,
            scalar=NULL, NBline=FALSE, nbins=100, log.axes="xy", xlab=NULL,
            ylab=NULL,  \dots)
binMeanVar(x, group, nbins=100, common.dispersion=FALSE, object=NULL)
}

\arguments{ 
\item{object}{\code{DGEList} object containing the raw data and dispersion value. According the method desired for computing the dispersion, either \code{estimateCommonDisp} and (possibly) \code{estimateTagwiseDisp} should be run on the \code{DGEList} object before using \code{plotMeanVar}. The argument \code{object} must be supplied in the function \code{binMeanVar} if common dispersion values are to be computed for each bin.}

\item{meanvar}{list (optional) containing the output from \code{binMeanVar} or the returned value of \code{plotMeanVar}. Providing this object as an argument will save time in computing the gene means and variances when producing a mean-variance plot. }

\item{show.raw.vars}{logical, whether or not to display the raw (pooled) genewise variances on the mean-variance plot. Default is \code{FALSE}.}

\item{show.tagwise.vars}{logical, whether or not to display the estimated genewise variances on the mean-variance plot (note that `tag' and `gene' are synonymous). Default is \code{FALSE}.}

\item{show.binned.common.disp.vars}{logical, whether or not to compute the common dispersion for each bin of genes and show the variances computed from those binned common dispersions and the mean expression level of the respective bin of genes. Default is \code{FALSE}.}

\item{show.ave.raw.vars}{logical, whether or not to show the average of the raw variances for each bin of genes plotted against the average expression level of the genes in the bin. Averages are taken on the square root scale as regular arithmetic means are likely to be upwardly biased for count data, whereas averaging on the square scale gives a better summary of the mean-variance relationship in the data. The default is \code{TRUE}.}

\item{scalar}{vector (optional) of scaling values to divide counts by. Would expect to have this the same length as the number of columns in the count matrix (i.e. the number of libraries).}

\item{NBline}{logical, whether or not to add a line on the graph showing the mean-variance relationship for a NB model with common dispersion.}

\item{nbins}{scalar giving the number of bins (formed by using the quantiles of the genewise mean expression levels) for  which to compute average means and variances for exploring the mean-variance relationship. Default is \code{100} bins}

\item{log.axes}{character vector indicating if any of the axes should use a log scale. Default is \code{"xy"}, which makes both y and x axes on the log scale. Other valid options are \code{"x"} (log scale on x-axis only), \code{"y"} (log scale on y-axis only) and \code{""} (linear scale on x- and y-axis).}

\item{xlab}{character string giving the label for the x-axis. Standard graphical parameter. If left as the default \code{NULL}, then the x-axis label will be set to "logConc".}

\item{ylab}{character string giving the label for the y-axis. Standard graphical parameter. If left as the default \code{NULL}, then the x-axis label will be set to "logConc".}

\item{\dots}{further arguments passed on to \code{plot}}

\item{x}{matrix of count data, with rows representing genes and columns representing samples}

\item{group}{factor giving the experimental group or condition to which each sample (i.e. column of \code{x} or element of {y}) belongs}

\item{common.dispersion}{logical, whether or not to compute the common dispersion for each bin of genes.}

}

\value{ \code{plotMeanVar} produces a mean-variance plot for the DGE data using the options described above. \code{plotMeanVar} and \code{binMeanVar} both return a list with the following components:
	\item{avemeans}{vector of the average expression level within each bin of genes, with the average taken on the square-root scale}
	\item{avevars}{vector of the average raw pooled gene-wise variance within each bin of genes, with the average taken on the square-root scale}
	\item{bin.means}{list containing the average (mean) expression level for genes divided into bins based on amount of expression}
	\item{bin.vars}{list containing the pooled variance for genes divided into bins based on amount of expression}
	\item{means}{vector giving the mean expression level for each gene}
	\item{vars}{vector giving the pooled variance for each gene}
	\item{bins}{list giving the indices of the genes in each bin, ordered from lowest expression bin to highest}
}

\details{
This function is useful for exploring the mean-variance relationship in the data. Raw variances are, for each gene, the pooled variance of the counts from each sample, divided by a scaling factor (by default the effective library size). The function will plot the average raw variance for genes split into \code{nbins} bins by overall expression level. The averages are taken on the square-root scale as for count data the arithmetic mean is upwardly biased. Taking averages on the square-root scale provides a useful summary of how the variance of the gene counts change with respect to expression level (abundance). A line showing the Poisson mean-variance relationship (mean equals variance) is always shown to illustrate how the genewise variances may differ from a Poisson mean-variance relationship. Optionally, the raw variances and estimated genewise variances can also be plotted. Estimated genewise variances can be calculated using either qCML estimates of the genewise dispersions (\code{estimateTagwiseDisp}) or Cox-Reid conditional inference estimates (\code{CRDisp}). A log-log scale is used for the plot.
}



\author{Davis McCarthy}
\examples{
y <- matrix(rnbinom(1000,mu=10,size=2),ncol=4)
d <- DGEList(counts=y,group=c(1,1,2,2),lib.size=c(1000:1003))
plotMeanVar(d) # Produce a straight-forward mean-variance plot
# Produce a mean-variance plot with the raw variances shown and save the means
# and variances for later use
meanvar <- plotMeanVar(d, show.raw.vars=TRUE) 
## If we want to show estimated genewise variances on the plot, we must first estimate them!
d <- estimateCommonDisp(d) # Obtain an estimate of the dispersion parameter
d <- estimateTagwiseDisp(d)  # Obtain genewise dispersion estimates
# Use previously saved object to speed up plotting
plotMeanVar(d, meanvar=meanvar, show.tagwise.vars=TRUE, NBline=TRUE) 
## We could also estimate common/genewise dispersions using the Cox-Reid methods with an
## appropriate design matrix
}

\seealso{
\code{\link{plotMeanVar2}} is related in purpose but works from raw counts and uses standardized residuals.

\code{\link{plotMDS.DGEList}}, \code{\link{plotSmear}} and
\code{\link{maPlot}} provide more ways of visualizing DGE data.

}

\concept{Dispersion estimation}
