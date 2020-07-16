\name{glmFit}
\alias{glmFit}
\alias{glmFit.DGEList}
\alias{glmFit.SummarizedExperiment}
\alias{glmFit.default}
\alias{glmLRT}

\title{Genewise Negative Binomial Generalized Linear Models}

\description{Fit a negative binomial generalized log-linear model to the read counts for each gene.
Conduct genewise statistical tests for a given coefficient or coefficient contrast.}

\usage{
\method{glmFit}{DGEList}(y, design=NULL, dispersion=NULL, prior.count=0.125, start=NULL, \dots)
\method{glmFit}{SummarizedExperiment}(y, design=NULL, dispersion=NULL, prior.count=0.125, start=NULL, \dots)
\method{glmFit}{default}(y, design=NULL, dispersion=NULL, offset=NULL, lib.size=NULL, weights=NULL,
       prior.count=0.125, start=NULL, \dots)
glmLRT(glmfit, coef=ncol(glmfit$design), contrast=NULL)
}

\arguments{
\item{y}{an object that contains the raw counts for each library (the measure of expression level); alternatively, a matrix of counts, or a \code{DGEList} object with (at least) elements \code{counts} (table of unadjusted counts) and \code{samples} (data frame containing information about experimental group, library size and normalization factor for the library size), or a \code{SummarizedExperiment} object with (at least) an component \code{counts} in its \code{assays}.}

\item{design}{numeric matrix giving the design matrix for the genewise linear models.
Must be of full column rank.
Defaults to a single column of ones, equivalent to treating the columns as replicate libraries.}

\item{dispersion}{numeric scalar, vector or matrix of negative binomial dispersions. Can be a common value for all genes, a vector of dispersion values with one for each gene, or a matrix of dispersion values with one for each observation. If \code{NULL} will be extracted from \code{y}, with order of precedence: genewise dispersion, trended dispersions, common dispersion.}

\item{offset}{numeric matrix of same size as \code{y} giving offsets for the log-linear models.  Can be a scalar or a vector of length \code{ncol(y)}, in which case it is expanded out to a matrix.}

\item{weights}{optional numeric matrix giving prior weights for the observations (for each library and gene) to be used in the GLM calculations.}

\item{lib.size}{numeric vector of length \code{ncol(y)} giving library sizes. Only used if \code{offset=NULL}, in which case \code{offset} is set to \code{log(lib.size)}. Defaults to \code{colSums(y)}.}

\item{prior.count}{average prior count to be added to observation to shrink the estimated log-fold-changes towards zero.}

\item{start}{optional numeric matrix of initial estimates for the linear model coefficients.}

\item{\dots}{other arguments are passed to lower level fitting functions.}

\item{glmfit}{a \code{DGEGLM} object, usually output from \code{glmFit}.}

\item{coef}{integer or character vector indicating which coefficients of the linear model are to be tested equal to zero.  Values must be columns or column names of \code{design}. Defaults to the last coefficient.  Ignored if \code{contrast} is specified.}

\item{contrast}{numeric vector or matrix specifying one or more contrasts of the linear model coefficients to be tested equal to zero. Number of rows must equal to the number of columns of \code{design}. If specified, then takes precedence over \code{coef}.}
}

\value{
\code{glmFit} produces an object of class \code{DGEGLM} containing components \code{counts}, \code{samples}, \code{genes} and \code{abundance} from \code{y} plus the following new components:
	\item{design}{design matrix as input.}
	\item{weights}{matrix of weights as input.}
	\item{df.residual}{numeric vector of residual degrees of freedom, one for each gene.}
	\item{offset}{numeric matrix of linear model offsets.}
	\item{dispersion}{vector of dispersions used for the fit.}
	\item{coefficients}{numeric matrix of estimated coefficients from the glm fits, on the natural log scale, of size \code{nrow(y)} by \code{ncol(design)}.}
	\item{unshrunk.coefficients}{numeric matrix of estimated coefficients from the glm fits when no log-fold-changes shrinkage is applied, on the natural log scale, of size \code{nrow(y)} by \code{ncol(design)}. It exists only when \code{prior.count} is not 0.}
	\item{fitted.values}{matrix of fitted values from glm fits, same number of rows and columns as \code{y}.}
	\item{deviance}{numeric vector of deviances, one for each gene.}

\code{glmLRT} produces objects of class \code{DGELRT} with the same components as for \code{glmFit} plus the following:
	\item{table}{data frame with the same rows as \code{y} containing the log2-fold-changes, likelhood ratio statistics and p-values, ready to be displayed by \code{topTags}.}
	\item{comparison}{character string describing the coefficient or the contrast being tested.}

The data frame \code{table} contains the following columns:
	\item{logFC}{log2-fold change of expression between conditions being tested.}
	\item{logCPM}{average log2-counts per million, the average taken over all libraries in \code{y}.}
	\item{LR}{likelihood ratio statistics.}
	\item{PValue}{p-values.}
}

\details{
\code{glmFit} and \code{glmLRT} implement generalized linear model (glm) methods developed by McCarthy et al (2012).

\code{glmFit} fits genewise negative binomial glms, all with the same design matrix but possibly different dispersions, offsets and weights.
When the design matrix defines a one-way layout, or can be re-parametrized to a one-way layout, the glms are fitting very quickly using \code{\link{mglmOneGroup}}.
Otherwise the default fitting method, implemented in \code{\link{mglmLevenberg}}, uses a Fisher scoring algorithm with Levenberg-style damping.

Positive \code{prior.count} cause the returned coefficients to be shrunk in such a way that fold-changes between the treatment conditions are decreased.
In particular, infinite fold-changes are avoided.
Larger values cause more shrinkage.
The returned coefficients are affected but not the likelihood ratio tests or p-values.

\code{glmLRT} conducts likelihood ratio tests for one or more coefficients in the linear model.
If \code{coef} is used, the null hypothesis is that all the coefficients indicated by \code{coef} are equal to zero.
If \code{contrast} is non-null, then the null hypothesis is that the specified contrasts of the coefficients are equal to zero.
For example, a contrast of \code{c(0,1,-1)}, assuming there are three coefficients, would test the hypothesis that the second and third coefficients are equal.
}

\references{
McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
\emph{Nucleic Acids Research} 40, 4288-4297.
\url{https://doi.org/10.1093/nar/gks042}
}

\author{Davis McCarthy and Gordon Smyth}

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

# Likelihood ratio tests for trend
results <- glmLRT(fit, coef=2)
topTags(results)

# Estimate the dispersion (may be unreliable with so few genes)
d <- estimateGLMCommonDisp(d, design, verbose=TRUE)
}

\seealso{
Low-level computations are done by \code{\link{mglmOneGroup}} or \code{\link{mglmLevenberg}}.

\code{\link{topTags}} displays results from \code{glmLRT}.
}

\concept{Model fit}
\concept{Differential expression}
