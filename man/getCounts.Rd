\name{getCounts}
\alias{getCounts}
\alias{getOffset}
\alias{getDispersion}

\title{Extract Specified Component of a DGEList Object}

\usage{
getCounts(y) 
getOffset(y) 
getDispersion(y)
}

\arguments{ 
\item{y}{\code{DGEList} object containing (at least) the elements \code{counts} (table of raw counts), \code{group} (factor indicating group) and \code{lib.size} (numeric vector of library sizes)}
}

\description{
\code{getCounts(y)} returns the matrix of read counts \code{y$counts}.

\code{getOffset(y)} returns offsets for the log-linear predictor account for sequencing depth and possibly other normalization factors.
Specifically it returns the matrix \code{y$offset} if it is non-null, otherwise it returns the log product of \code{lib.size} and \code{norm.factors} from \code{y$samples}.

\code{getDispersion(y)} returns the most complex dispersion estimates (common, trended or genewise) found in \code{y}.
}

\value{\code{getCounts} returns the matrix of counts.
\code{getOffset} returns a numeric matrix or vector.
\code{getDispersion} returns vector of dispersion values.
}

\author{Mark Robinson, Davis McCarthy, Gordon Smyth}

\examples{
# generate raw counts from NB, create list object
y <- matrix(rnbinom(20,size=5,mu=10),5,4)
d <- DGEList(counts=y, group=c(1,1,2,2), lib.size=1001:1004)
getCounts(d)
getOffset(d)
d <- estimateCommonDisp(d)
getDispersion(d)
}

\seealso{\code{\link[edgeR:DGEList-class]{DGEList-class}}}

\concept{edgeR classes}
