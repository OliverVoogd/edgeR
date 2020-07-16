\name{condLogLikDerSize}
\alias{condLogLikDerSize}
\alias{condLogLikDerDelta}

\title{Conditional Log-Likelihood of the Dispersion for a Single Group of Replicate Libraries}

\description{Derivatives of the negative-binomial log-likelihood with respect to the dispersion parameter for each gene, conditional on the mean count, for a single group of replicate libraries of the same size.}

\usage{
condLogLikDerSize(y, r, der=1L)
condLogLikDerDelta(y, delta, der=1L)
}

\arguments{ 
\item{y}{matrix of counts, all counts in each row having the same population mean}
\item{r}{numeric vector or scalar, size parameter of negative binomial distribution, equal to 1/dispersion}
\item{delta}{numeric vector or scalar, delta parameter of negative binomial, equal to dispersion/(1+dispersion)}
\item{der}{integer specifying derivative required, either 0 (the function), 1 (first derivative) or 2 (second derivative)}
}

\value{vector of row-wise derivatives with respect to \code{r} or \code{delta}}

\details{The library sizes must be equalized before running this function. This function carries out the actual mathematical computations for the conditional log-likelihood and its derivatives, calculating the conditional log-likelihood for each gene.
Derivatives are with respect to either the size (\code{r}) or the delta parametrization (\code{delta}) of the dispersion.
}


\author{Mark Robinson, Davis McCarthy, Gordon Smyth}

\examples{
y <- matrix(rnbinom(10,size=1,mu=10),nrow=5)
condLogLikDerSize(y,r=1,der=1)
condLogLikDerDelta(y,delta=0.5,der=1)
}

\concept{Dispersion estimation}
