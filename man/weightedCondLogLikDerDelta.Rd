\name{weightedCondLogLikDerDelta}
\alias{weightedCondLogLikDerDelta}

\title{Weighted Conditional Log-Likelihood in Terms of Delta}

\description{Weighted conditional log-likelihood parameterized in terms of delta (\code{phi / (phi+1)}) for a given gene, maximized to find the smoothed (moderated) estimate of the dispersion parameter}

\usage{
weightedCondLogLikDerDelta(y, delta, tag, prior.n=10, ntags=nrow(y[[1]]), der=0)
}

\arguments{
\item{y}{list with elements comprising the matrices of count data (or pseudocounts) for the different groups}

\item{delta}{delta (\code{phi / (phi+1)})parameter of negative binomial}

\item{tag}{gene at which the weighted conditional log-likelihood is evaluated}

\item{prior.n}{smoothing paramter that indicates the weight to put on the common likelihood compared to the individual gene's likelihood; default \code{10} means that the common likelihood is given 10 times the weight of the individual gene's likelihood in the estimation of the genewise dispersion}

\item{ntags}{numeric scalar number of genes in the dataset to be analysed}

\item{der}{derivative, either 0 (the function), 1 (first derivative) or 2 (second derivative)}
}

\value{ numeric scalar of function/derivative evaluated for the  given gene and delta}

\details{
This function computes the weighted conditional log-likelihood for a given gene, parameterized in terms of delta. The value of delta that maximizes the weighted conditional log-likelihood is converted back to the \code{phi} scale, and this value is the estimate of the smoothed (moderated) dispersion parameter for that particular gene. The delta scale for convenience (delta is bounded between 0 and 1). 
Users should note that `tag' and `gene' are synonymous when interpreting the names of the arguments for this function.
}


\author{Mark Robinson, Davis McCarthy}

\examples{
counts<-matrix(rnbinom(20,size=1,mu=10),nrow=5)
d<-DGEList(counts=counts,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
y<-splitIntoGroups(d)
ll1<-weightedCondLogLikDerDelta(y,delta=0.5,tag=1,prior.n=10,der=0)
ll2<-weightedCondLogLikDerDelta(y,delta=0.5,tag=1,prior.n=10,der=1)
}

\concept{Dispersion estimation}
