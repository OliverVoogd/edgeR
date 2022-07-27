\name{getPriorN}
\alias{getPriorN}

\title{Get a Recommended Value for Prior N from DGEList Object}

\description{Returns the \code{lib.size} component of the \code{samples} component of  DGEList object multiplied by the \code{norm.factors} component}

\usage{
getPriorN(y, design=NULL, prior.df=20)
}
\arguments{ 

\item{y}{a \code{DGEList} object with (at least) elements \code{counts} (table of unadjusted counts) and \code{samples} (data frame containing information about experimental group, library size and normalization factor for the library size)}

\item{design}{numeric matrix (optional argument) giving the design matrix for the GLM that is to be fit. Must be of full column rank. If provided \code{design} is used to determine the number of parameters to be fit in the statistical model and therefore the residual degrees of freedom. If left as the default (\code{NULL}) then the \code{y$samples$group} element of the \code{DGEList} object is used to determine the residual degrees of freedom.}

\item{prior.df}{numeric scalar giving the weight, in terms of prior degrees of freedom, to be given to the common parameter likelihood when estimating genewise dispersion estimates.}

}

\value{\code{getPriorN} returns a numeric scalar }

\details{
When estimating genewise dispersion values using \code{\link{estimateTagwiseDisp}} or \code{\link{estimateGLMTagwiseDisp}} we need to decide how much weight to give to the common parameter likelihood in order to smooth (or stabilize) the dispersion estimates. The best choice of value for the \code{prior.n} parameter varies between datasets depending on the number of samples in the dataset and the complexity of the model to be fit. The value of \code{prior.n} should be inversely proportional to the residual degrees of freedom. We have found that choosing a value for \code{prior.n} that is equivalent to giving the common parameter likelihood 20 degrees of freedom generally gives a good amount of smoothing for the genewise dispersion estimates. This function simply recommends an appropriate value for \code{prior.n}---to be used as an argument for \code{\link{estimateTagwiseDisp}} or \code{\link{estimateGLMTagwiseDisp}}---given the experimental design at hand and the chosen prior degrees of freedom.
}

\author{Davis McCarthy, Gordon Smyth}

\examples{
# generate raw counts from NB, create list object
y<-matrix(rnbinom(20,size=1,mu=10),nrow=5)
d<-DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
getPriorN(d)
}

\seealso{
\code{\link{DGEList}} for more information about the \code{DGEList} class.
\code{\link{as.matrix.DGEList}}.
}

\concept{Dispersion estimation}
