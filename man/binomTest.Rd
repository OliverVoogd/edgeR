\name{binomTest}
\alias{binomTest}

\title{Exact Binomial Tests for Comparing Two Digital Libraries}
\description{
Computes p-values for differential abundance for each gene between two digital libraries,
conditioning on the total count for each gene.
The counts in each group as a proportion of the whole are assumed to follow a binomial distribution.
}

\usage{
binomTest(y1, y2, n1=sum(y1), n2=sum(y2), p=n1/(n1+n2))
}

\arguments{
\item{y1}{integer vector giving the count for each gene in the first library.
Non-integer values are rounded to the nearest integer.}
\item{y2}{integer vector giving the count for each gene in the second library.
Of same length as \code{y1}.
Non-integer values are rounded to the nearest integer.}
\item{n1}{total number of counts in the first library, across all genes.
Non-integer values are rounded to the nearest integer. Not required if \code{p} is supplied.}
\item{n2}{total number of counts in the second library, across all genes.
Non-integer values are rounded to the nearest integer. Not required if \code{p} is supplied.}
\item{p}{expected proportion of \code{y1} to the total for each gene under the null hypothesis.}
}

\details{
This function can be used to compare two libraries from SAGE, RNA-Seq, ChIP-Seq or other sequencing technologies with respect to technical variation.

An exact two-sided binomial test is computed for each gene.
This test is closely related to Fisher's exact test for 2x2 contingency tables but, unlike Fisher's test, it conditions on the total number of counts for each gene.
The null hypothesis is that the expected counts are in the same proportions as the library sizes, i.e., that the binomial probability for the first library is \code{n1/(n1+n2)}.

The two-sided rejection region is chosen analogously to Fisher's test.
Specifically, the rejection region consists of those values with smallest probabilities
under the null hypothesis.

When the counts are reasonably large, the binomial test, Fisher's test and Pearson's chisquare all give the same results.
When the counts are smaller, the binomial test is usually to be preferred in this context.

This function replaces the earlier \code{sage.test} functions in the statmod and sagenhaft packages.
It produces the same results as \code{\link{binom.test}} in the stats packge, but is much faster.
}

\value{
Numeric vector of p-values.
}

\references{
\url{http://en.wikipedia.org/wiki/Binomial_test}

\url{http://en.wikipedia.org/wiki/Fisher's_exact_test}

\url{http://en.wikipedia.org/wiki/Serial_analysis_of_gene_expression}

http://en.wikipedia.org/wiki/RNA-Seq
}

\author{Gordon Smyth}

\seealso{
\code{\link[statmod:sage.test]{sage.test}} (statmod package), \code{\link{binom.test}} (stats package)
}

\examples{
binomTest(c(0,5,10),c(0,30,50),n1=10000,n2=15000)
#  Univariate equivalents:
binom.test(5,5+30,p=10000/(10000+15000))$p.value
binom.test(10,10+50,p=10000/(10000+15000))$p.value
}

\concept{Differential expression}
