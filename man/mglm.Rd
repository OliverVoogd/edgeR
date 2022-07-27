\name{mglm}
\alias{mglm}
\alias{mglmOneGroup}
\alias{mglmOneWay}
\alias{mglmLevenberg}
\alias{designAsFactor}

\title{Fit Negative Binomial Generalized Linear Models to Multiple Response Vectors: Low Level Functions}

\description{
Fit the same log-link negative binomial or Poisson generalized linear model (GLM) to each row of a matrix of counts.
}

\usage{
mglmOneGroup(y, dispersion = 0, offset = 0, weights = NULL,
             coef.start = NULL, maxit = 50, tol = 1e-10, verbose = FALSE)
mglmOneWay(y, design = NULL, group = NULL, dispersion = 0, offset = 0, weights = NULL,
             coef.start = NULL, maxit = 50, tol = 1e-10)
mglmLevenberg(y, design, dispersion = 0, offset = 0, weights = NULL,
             coef.start = NULL, start.method = "null", maxit = 200, tol = 1e-06)
designAsFactor(design)
}

\arguments{
\item{y}{numeric matrix containing the negative binomial counts.  Rows for genes and columns for libraries.}

\item{design}{numeric matrix giving the design matrix of the GLM.
Assumed to be full column rank.
This is a required argument for \code{mglmLevernberg} and \code{code{designAsFactor}}.
For \code{mglmOneWay}, it defaults to \code{model.matrix(~0+group)} if \code{group} is specified and otherwise to \code{model.matrix( ~1)}.
}

\item{group}{factor giving group membership for oneway layout.
If both \code{design} and \code{group} are both specified, then they must agree in terms of \code{designAsFactor}.
If \code{design=NULL}, then a group-means design matrix is implied.}

\item{dispersion}{numeric scalar or vector giving the dispersion parameter for each GLM.
Can be a scalar giving one value for all genes, or a vector of length equal to the number of genes giving genewise dispersions.}

\item{offset}{numeric vector or matrix giving the offset that is to be included in the log-linear model predictor.  Can be a scalar, a vector of length equal to the number of libraries, or a matrix of the same size as \code{y}.}

\item{weights}{numeric vector or matrix of non-negative quantitative weights.
Can be a vector of length equal to the number of libraries, or a matrix of the same size as \code{y}.}

\item{coef.start}{numeric matrix of starting values for the linear model coefficients.
Number of rows should agree with \code{y} and number of columns should agree with \code{design}.
For \code{mglmOneGroup}, a numeric vector or a matrix with one column.
This argument does not usually need to be set as the automatic starting values perform well.}

\item{start.method}{method used to generate starting values when \code{coef.stat=NULL}. Possible values are \code{"null"} to start from the null model of equal expression levels or \code{"y"} to use the data as starting value for the mean.}

\item{tol}{numeric scalar giving the convergence tolerance. For \code{mglmOneGroup}, convergence is judged successful when the step size falls below \code{tol} in absolute size.}

\item{maxit}{integer giving the maximum number of iterations for the Fisher scoring algorithm. The iteration will be stopped when this limit is reached even if the convergence criterion hasn't been satisfied.}

\item{verbose}{logical. If \code{TRUE}, warnings will be issued when \code{maxit} iterations are exceeded before convergence is achieved.}
}

\details{
These functions are low-level work-horses used by higher-level functions in the edgeR package, especially by \code{glmFit}.
Most users will not need to call these functions directly.

The functions \code{mglmOneGroup}, \code{mglmOneWay} and \code{mglmLevenberg} all fit a negative binomial GLM to each row of \code{y}.
The row-wise GLMS all have the same design matrix but possibly different dispersions, offsets and weights.
These functions are all low-level in that they operate on atomic objects (numeric matrices and vectors).

\code{mglmOneGroup} fits an intercept only model to each response vector.
In other words, it treats all the libraries as belonging to one group.
It implements Fisher scoring with a score-statistic stopping criterion for each gene.
Excellent starting values are available for the null model so this function seldom has any problems with convergence.
It is used by other edgeR functions to compute the overall abundance for each gene.

\code{mglmOneWay} fits a oneway layout to each response vector.
It treats the libraries as belonging to a number of groups and calls \code{mglmOneGroup} for each group.

\code{mglmLevenberg} fits an arbitrary log-linear model to each response vector.
It implements a Levenberg-Marquardt modification of the GLM scoring algorithm to prevent divergence.
The main computation is implemented in C++.

All these functions treat the dispersion parameter of the negative binomial distribution as a known input.

\code{designAsFactor} is used to convert a general design matrix into a oneway layout if that is possible.
It determines how many distinct row values the design matrix is capable of computing and returns a factor with a level for each possible distinct value.
}

\value{
\code{mglmOneGroup} produces a numeric vector of coefficients.

\code{mglmOneWay} produces a list with the following components: 
	\item{coefficients}{matrix of estimated coefficients for the linear models. Rows correpond to rows of \code{y} and columns to columns of \code{design}.}
	\item{fitted.values}{matrix of fitted values. Of same dimensions as \code{y}.}

\code{mglmLevenberg} produces a list with the following components: 
	\item{coefficients}{matrix of estimated coefficients for the linear models.}
	\item{fitted.values}{matrix of fitted values.}
	\item{deviance}{numeric vector of residual deviances.}
	\item{iter}{number of iterations used.}
	\item{fail}{logical vector indicating genes for which the maximum damping was exceeded before convergence was achieved.}

\code{designAsFactor} returns a factor of length equal to \code{nrow(design)}.
}

\references{
McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
\emph{Nucleic Acids Research} 40, 4288-4297.
\doi{10.1093/nar/gks042}
}

\author{Gordon Smyth, Yunshun Chen, Davis McCarthy, Aaron Lun.  C++ code by Aaron Lun.}

\seealso{
Most users will call either \code{\link{glmFit}}, the higher-level function offering more object-orientated GLM modelling of DGE data, or else \code{\link{exactTest}}, which is designed for oneway layouts.
}

\examples{
y <- matrix(rnbinom(1000, mu = 10, size = 2), ncol  =  4)
lib.size <- colSums(y)
dispersion <- 0.1

## Compute intercept for each row
beta <- mglmOneGroup(y, dispersion = dispersion, offset = log(lib.size))

## Unlogged intercepts add to one:
sum(exp(beta))

## Fit the NB GLM to the counts with a given design matrix
f1 <- factor(c(1,1,2,2))
f2 <- factor(c(1,2,1,2))
X <- model.matrix(~ f1 + f2)
fit <- mglmLevenberg(y, X, dispersion = dispersion, offset = log(lib.size))
head(fit$coefficients)
}

\concept{Model fit}
