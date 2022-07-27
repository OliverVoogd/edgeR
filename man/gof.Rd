\name{gof}
\alias{gof}

\title{Goodness of Fit Tests for Multiple GLM Fits}

\description{Conducts deviance goodness of fit tests for each fit in a \code{DGEGLM} object}

\usage{
gof(glmfit, pcutoff = 0.1, adjust = "holm", plot = FALSE,
    main = "qq-plot of residual deviances", \dots)
}

\arguments{ 
\item{glmfit}{a \code{DGEGLM} object containing results from fitting NB GLMs to genes in a DGE dataset with a global dispersion model. Usually this is output from \code{glmFit}.}
\item{pcutoff}{scalar giving the cut-off value for the Holm-adjusted p-value. Genes with Holm-adjusted p-values lower than this cutoff value are flagged as `dispersion outlier' genes.}
\item{adjust}{method used to adjust goodness of fit p-values for multiple testing.}
\item{plot}{logical, if \code{TRUE} a qq-plot is produced.}
\item{main}{character, title for the plot.}
\item{\dots}{other arguments are passed to \code{qqnorm}.}
}

\details{
This function is useful for evaluating the adequacy of a global dispersion model, such as a constant or trended dispersion.
If \code{plot=TRUE}, then it produces a qq-plot similar to those in Figure 2 of McCarthy et al (2012).
}

\note{
This function should not be used with tagwise estimated dispersions such as those from \code{\link{estimateGLMTagwiseDisp}} or \code{\link{estimateDisp}}.
\code{glmfit} should contain trended or constant dispersions.
}

\value{
A list with the following components:
\item{gof.statistics}{numeric vector of deviance statistics, which are the statistics used for the goodness of fit test}
\item{gof.pvalues}{numeric vector of p-values providing evidence of poor fit; computed from the chi-square distribution on the residual degrees of freedom from the GLM fits.}
\item{outlier}{logical vector indicating whether or not each gene is a `dispersion outlier' (i.e., the model fit is poor for that gene indicating that the dispersion estimate is not good for that gene).} 
\item{df}{scalar, the residual degrees of freedom from the GLM fit for which the goodness of fit statistics have been computed. Also the degrees of freedom for the goodness of fit statistics for the LR (chi-quare) test for significance.}

If \code{plot=TRUE}, then a plot is also produced on the current graphics device.
}

\author{Davis McCarthy and Gordon Smyth}

\references{
McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
\emph{Nucleic Acids Research} 40, 4288-4297
\doi{10.1093/nar/gks042}
}

\examples{
nlibs <- 3
ngenes <- 100
dispersion.true <- 0.1

# Make first gene respond to covariate x
x <- 0:2
design <- model.matrix(~x)
beta.true <- cbind(Beta1=2,Beta2=c(2,rep(0,ngenes-1)))
mu.true <- 2^(beta.true \%*\% t(design))

# Generate count data
y <- rnbinom(ngenes*nlibs,mu=mu.true,size=1/dispersion.true)
y <- matrix(y,ngenes,nlibs)
colnames(y) <- c("x0","x1","x2")
rownames(y) <- paste("gene",1:ngenes,sep=".")
d <- DGEList(y)

# Normalize
d <- calcNormFactors(d)

# Fit the NB GLMs
fit <- glmFit(d, design, dispersion=dispersion.true)
# Check how good the fit is for each gene
gof(fit)
}

\seealso{
\code{\link{qqnorm}}.

\code{\link{glmFit}} for more information on fitting NB GLMs to DGE data.
}

\concept{Dispersion estimation}
