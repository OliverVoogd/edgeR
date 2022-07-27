\name{commonCondLogLikDerDelta}
\alias{commonCondLogLikDerDelta}

\title{Conditional Log-Likelihoods in Terms of Delta}

\description{Common conditional log-likelihood parameterized in terms of delta (\code{phi / (phi+1)}) }

\usage{
commonCondLogLikDerDelta(y, delta, der = 0) 
}
\arguments{ 

\item{y}{list with elements comprising the matrices of count data (or pseudocounts) for the different groups}

\item{delta}{delta (\code{phi / (phi+1)}) parameter of negative binomial}

\item{der}{derivative, either 0 (the function), 1 (first derivative) or 2 (second derivative)}
}

\value{ numeric scalar of function/derivative evaluated at given delta}


\details{
The common conditional log-likelihood is constructed by summing over all of the individual genewise conditional log-likelihoods. The common conditional log-likelihood is taken as a function of the dispersion parameter (\code{phi}), and here parameterized in terms of delta (\code{phi / (phi+1)}). The value of delta that maximizes the common conditional log-likelihood is converted back to the \code{phi} scale, and this value is the estimate of the common dispersion parameter used by all genes.
}

\author{Davis McCarthy}

\examples{
counts<-matrix(rnbinom(20,size=1,mu=10),nrow=5)
d<-DGEList(counts=counts,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
y<-splitIntoGroups(d)
ll1<-commonCondLogLikDerDelta(y,delta=0.5,der=0)
ll2<-commonCondLogLikDerDelta(y,delta=0.5,der=1)
}

\seealso{
\code{\link{estimateCommonDisp}} is the user-level function for estimating the common dispersion parameter.
}

\concept{Dispersion estimation}
