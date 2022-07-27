\name{estimateGLMTagwiseDisp}
\alias{estimateGLMTagwiseDisp}
\alias{estimateGLMTagwiseDisp.DGEList}
\alias{estimateGLMTagwiseDisp.default}


\title{Empirical Bayes Tagwise Dispersions for Negative Binomial GLMs}

\description{
Compute an empirical Bayes estimate of the negative binomial dispersion parameter for each tag, with expression levels specified by a log-linear model.
}

\usage{
\S3method{estimateGLMTagwiseDisp}{DGEList}(y, design=NULL, prior.df=10,
            trend=!is.null(y$trended.dispersion), span=NULL, \dots)
\S3method{estimateGLMTagwiseDisp}{default}(y, design=NULL, offset=NULL, dispersion,
            prior.df=10, trend=TRUE, span=NULL, AveLogCPM=NULL,
            weights=NULL, \dots)
}

\arguments{
\item{y}{matrix of counts or a \code{DGEList} object.}
\item{design}{numeric design matrix, as for \code{\link{glmFit}}.}
\item{trend}{logical. Should the prior be the trended dispersion (\code{TRUE}) or the common dispersion (\code{FALSE})?}
\item{offset}{offset matrix for the log-linear model, as for \code{\link{glmFit}}.  Defaults to the log-effective library sizes.}
\item{dispersion}{common or trended dispersion estimates, used as an initial estimate for the tagwise estimates.}
\item{prior.df}{prior degrees of freedom.}
\item{span}{width of the smoothing window, in terms of proportion of the data set. Default value decreases with the number of tags.}
\item{AveLogCPM}{numeric vector giving average log2 counts per million for each tag}
\item{weights}{optional numeric matrix giving observation weights}
\item{\ldots}{other arguments are passed to \code{\link{dispCoxReidInterpolateTagwise}}.}
}

\value{
\code{estimateGLMTagwiseDisp.DGEList} produces a \code{DGEList} object, which contains the tagwise dispersion parameter estimate for each tag for the negative binomial model that maximizes the Cox-Reid adjusted profile likelihood. The tagwise dispersions are simply added to the \code{DGEList} object provided as the argument to the function.

\code{estimateGLMTagwiseDisp.default} returns a vector of the tagwise dispersion estimates.
}

\details{
This function implements the empirical Bayes strategy proposed by McCarthy et al (2012)
for estimating the tagwise negative binomial dispersions.
The experimental conditions are specified by design matrix allowing for multiple explanatory factors.
The empirical Bayes posterior is implemented as a conditional likelihood with tag-specific weights,
and the conditional likelihood is computed using Cox-Reid approximate conditional likelihood (Cox and Reid, 1987).

The prior degrees of freedom determines the weight given to the global dispersion trend.
The larger the prior degrees of freedom, the more the tagwise dispersions are squeezed towards the global trend.

Note that the terms `tag' and `gene' are synonymous here. 
The function is only named `Tagwise' for historical reasons.

This function calls the lower-level function \code{\link{dispCoxReidInterpolateTagwise}}.
}

\references{
Cox, DR, and Reid, N (1987). Parameter orthogonality and approximate conditional inference.
\emph{Journal of the Royal Statistical Society Series B} 49, 1-39.

McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
\emph{Nucleic Acids Research} 40, 4288-4297.
\doi{10.1093/nar/gks042}
}

\author{Gordon Smyth, Davis McCarthy}
\examples{
y <- matrix(rnbinom(1000,mu=10,size=10),ncol=4)
d <- DGEList(counts=y,group=c(1,1,2,2),lib.size=c(1000:1003))
design <- model.matrix(~group, data=d$samples) # Define the design matrix for the full model
d <- estimateGLMTrendedDisp(d, design, min.n=10)
d <- estimateGLMTagwiseDisp(d, design)
summary(d$tagwise.dispersion)
}

\seealso{
\code{\link{estimateGLMCommonDisp}} for common dispersion or \code{\link{estimateGLMTrendedDisp}} for trended dispersion in the context of a generalized linear model.

\code{\link{estimateCommonDisp}} for common dispersion or \code{\link{estimateTagwiseDisp}} for tagwise dispersions in the context of a multiple group experiment (one-way layout).
}

\concept{Dispersion estimation}
