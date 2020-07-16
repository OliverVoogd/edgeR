\name{addPriorCount}
\alias{addPriorCount}

\title{Add a prior count}
\description{Add a library size-adjusted prior count to each observation.}

\usage{
addPriorCount(y, lib.size=NULL, offset=NULL, prior.count=1)
}

\arguments{
\item{y}{a numeric count matrix, with rows corresponding to genes and columns to libraries.}
\item{lib.size}{a numeric vector of library sizes.}
\item{offset}{a numeric vector or matrix of offsets.}
\item{prior.count}{a numeric scalar or vector of prior counts to be added to each gene.}
}

\details{
This function adds a positive prior count to each observation, often useful for avoiding zeroes during calculation of log-values.
For example, \code{\link{predFC}} will call this function to calculate shrunken log-fold changes.
\code{\link{aveLogCPM}} and \code{\link{cpm}} also use the same underlying code to calculate (average) log-counts per million.

The actual value added to the counts for each library is scaled according to the library size.
This ensures that the relative contribution of the prior is the same for each library.
Otherwise, a fixed prior would have little effect on a large library, but a big effect for a small library. 

The library sizes are also modified, with twice the scaled prior being added to the library size for each library.
To understand the motivation for this, consider that each observation is, effectively, a proportion of the total count in the library.
The addition scheme implemented here represents an empirical logistic transform and ensures that the proportion can never be zero or one.

If \code{offset} is supplied, this is used in favour of \code{lib.size} where \code{exp(offset)} is defined as the vector/matrix of library sizes.
If an offset matrix is supplied, this will lead to gene-specific scaling of the prior as described above.

Most use cases of this function will involve supplying a constant value to \code{prior.count} for all genes.
However, it is also possible to use gene-specific values by supplying a vector of length equal to the number of rows in \code{y}.
}

\value{
A list is returned containing \code{y}, a matrix of counts with the added priors; and \code{offset}, a CompressedMatrix containing the (log-transformed) modified library sizes.
}

\examples{
original <- matrix(rnbinom(1000, mu=20, size=10), nrow=200)
head(original)

out <- addPriorCount(original)
head(out$y)
head(out$offset)
}

\author{
Aaron Lun
}

\seealso{
\code{\link{aveLogCPM}},
\code{\link{cpm}},
\code{\link{predFC}}
}
