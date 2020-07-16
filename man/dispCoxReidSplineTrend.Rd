\name{dispCoxReidSplineTrend}
\alias{dispCoxReidSplineTrend}
\alias{dispCoxReidPowerTrend}

\title{Estimate Dispersion Trend for Negative Binomial GLMs}

\description{
Estimate trended dispersion parameters across multiple negative binomial generalized linear models using Cox-Reid adjusted profile likelihood.
}

\usage{
dispCoxReidSplineTrend(y, design, offset=NULL, df = 5, subset=10000, AveLogCPM=NULL,
                       method.optim="Nelder-Mead", trace=0)
dispCoxReidPowerTrend(y, design, offset=NULL, subset=10000, AveLogCPM=NULL,
                       method.optim="Nelder-Mead", trace=0)
}

\arguments{ 

\item{y}{numeric matrix of counts}

\item{design}{numeric matrix giving the design matrix for the GLM that is to be fit.}

\item{offset}{numeric scalar, vector or matrix giving the offset (in addition to the log of the effective library size) that is to be included in the NB GLM for the genes. If a scalar, then this value will be used as an offset for all genes and libraries. If a vector, it should be have length equal to the number of libraries, and the same vector of offsets will be used for each gene. If a matrix, then each library for each gene can have a unique offset, if desired. In \code{adjustedProfileLik} the \code{offset} must be a matrix with the same dimension as the table of counts.}

\item{df}{integer giving the degrees of freedom of the spline function, see \code{ns} in the splines package.}

\item{subset}{integer, number of rows to use in the calculation.  Rows used are chosen evenly spaced by AveLogCPM using \code{\link{cutWithMinN}}.}

\item{AveLogCPM}{numeric vector giving average log2 counts per million for each gene.}

\item{method.optim}{the method to be used in \code{optim}. See \code{\link{optim}} for more detail.}

\item{trace}{logical, should iteration information be output?}
}

\value{
List containing numeric vectors \code{dispersion} and \code{abundance} containing the estimated dispersion and abundance for each gene.
The vectors are of the same length as \code{nrow(y)}.
}

\details{
In the \code{edgeR} context, these are low-level functions called by \code{estimateGLMTrendedDisp}.

\code{dispCoxReidSplineTrend} and \code{dispCoxReidPowerTrend} fit abundance trends to the genewise dispersions.
\code{dispCoxReidSplineTrend} fits a regression spline whereas \code{dispCoxReidPowerTrend} fits a log-linear trend of the form \code{a*exp(abundance)^b+c}.
In either case, \code{optim} is used to maximize the adjusted profile likelihood (Cox and Reid, 1987).
}

\references{
Cox, DR, and Reid, N (1987). Parameter orthogonality and approximate conditional inference. \emph{Journal of the Royal Statistical Society Series B} 49, 1-39.
}

\author{Yunshun Chen, Davis McCarthy, Gordon Smyth}

\examples{
design <- matrix(1,4,1)
y <- matrix((rnbinom(400,mu=100,size=5)),100,4)
d1 <- dispCoxReidSplineTrend(y, design, df=3)
d2 <- dispCoxReidPowerTrend(y, design)
with(d2,plot(AveLogCPM,sqrt(dispersion)))
}

\seealso{
\code{\link{estimateGLMTrendedDisp}}
}

\concept{Dispersion estimation}
