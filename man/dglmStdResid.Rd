\name{plotMeanVar2}
\alias{plotMeanVar2}

\title{Plot Mean-Variance Relationship in DGE Data Using Standardized Residuals}

\description{
Group observations by size of the fitted value and plot average squared residual vs average fitted value.
}

\usage{
plotMeanVar2(y, design = NULL, dispersion = 0, offset = 0, nbins = 100, make.plot = TRUE,
          xlab = "Mean", ylab = "Ave. binned standardized residual", \dots)
}

\arguments{ 
  \item{y}{numeric matrix of counts, each row represents one genes, each column represents one DGE library.}

  \item{design}{numeric matrix giving the design matrix of the GLM. Assumed to be full column rank. Defaults to a intercept column.}

  \item{dispersion}{numeric scalar or vector giving the dispersion parameter for each GLM. Can be a scalar giving one value for all genes, or a vector of length equal to the number of genes giving genewise dispersions.}

  \item{offset}{numeric vector or matrix giving the offset that is to be included in teh log-linear model predictor. Can be a vector of length equal to the number of libraries, or a matrix of the same size as \code{y}.}

  \item{nbins}{scalar giving the number of bins (formed by using the quantiles of the genewise mean expression levels) for  which to compute average means and variances for exploring the mean-variance relationship. Default is \code{100} bins}

  \item{make.plot}{logical, whether or not to plot the mean standardized residual for binned data (binned on expression level). Provides a visualization of the mean-variance relationship. Default is \code{TRUE}.}

  \item{xlab}{character string giving the label for the x-axis. Standard graphical parameter. If left as the default, then the x-axis label will be set to "Mean".}

  \item{ylab}{character string giving the label for the y-axis. Standard graphical parameter. If left as the default, then the y-axis label will be set to "Ave. binned standardized residual".}

  \item{\dots}{other arguments are passed to \code{plot}}
}

\value{
Produces a mean-variance plot and returns a list with the following component:
	\item{mean}{numeric vector of average fitted values for bins}
	\item{var}{numeric vector of the average squared residuals for bins}
}

\details{
This function explores the mean-variance relationship in count data.
The function fits a Poisson or NB GLM model to each gene using the appropriate design matrix and computes squared ordinary residuals.
The residuals are standardized by the leverages but not by model variances.
The fitted values are divided into \code{nbins} bins and the mean fitted value and mean squared residual is computed for each bin.

This function is similar in purpose to \code{plotMeanVar} but uses standardized residuals instead of pooled variances.
}

\author{Davis McCarthy and Gordon Smyth}
\examples{
# Example with Poisson data
log2mu <- seq(from=0,to=10,length=1000)
y <- matrix(rpois(4*1000,lambda=2^log2mu),1000,4)
binned <- plotMeanVar2(y)
abline(0,1)
}

\seealso{
\code{\link{plotMeanVar}}.

\code{\link{plotMDS.DGEList}}, \code{\link{plotSmear}}, \code{\link{plotMD.DGEList}} and \code{\link{plotBCV}} provide other ways to visualize DGE data.
}

\concept{Dispersion estimation}
