\name{glmTreat}
\alias{glmTreat}

\title{Test for Differential Expression Relative to a Threshold}

\description{Conduct genewise statistical tests for a given coefficient or contrast relative to a specified fold-change threshold.}

\usage{
glmTreat(glmfit, coef = ncol(glmfit$design), contrast = NULL, lfc = log2(1.2),
         null = "interval")
}

\arguments{
\item{glmfit}{a \code{DGEGLM} object, usually output from \code{glmFit} or \code{glmQLFit}.}

\item{coef}{integer or character vector indicating which coefficients of the linear model are to be tested equal to zero.  Values must be columns or column names of \code{design}. Defaults to the last coefficient. Ignored if \code{contrast} is specified.}

\item{contrast}{numeric vector specifying the contrast of the linear model coefficients to be tested against the log2-fold-change threshold. Length must equal to the number of columns of \code{design}. If specified, then takes precedence over \code{coef}.}

\item{lfc}{numeric scalar specifying the absolute value of the log2-fold change threshold above which differential expression is to be considered.}

\item{null}{character string, choices are \code{"worst.case"} or \code{"interval"}.
If \code{"worst.case"}, then the null hypothesis asssumes that the true logFC is on the boundary of the possible values, either at \code{lfc} or \code{-lfc}, whichever gives the largest p-value.
This gives the most conservative results.
If \code{"interval"}, then the null hypotheses assumes the true logFC to belong to a bounded interval of possible values.}
}

\value{
\code{glmTreat} produces an object of class \code{DGELRT} with the same components as for \code{glmfit} plus the following:
	\item{lfc}{absolute value of the specified log2-fold-change threshold.}
	\item{table}{data frame with the same rows as \code{glmfit} containing the log2-fold-changes, average log2-counts per million and p-values, ready to be displayed by \code{topTags}.}
	\item{comparison}{character string describing the coefficient or the contrast being tested.}

The data frame \code{table} contains the following columns:
	\item{logFC}{shrunk log2-fold-change of expression between conditions being tested.}
	\item{unshrunk.logFC}{unshrunk log2-fold-change of expression between conditions being tested. Exists only when \code{prior.count} is not equal to 0 for \code{glmfit}.}
	\item{logCPM}{average log2-counts per million, the average taken over all libraries.}
	\item{PValue}{p-values.}
}

\details{
\code{glmTreat} implements a test for differential expression relative to a minimum required fold-change threshold.
Instead of testing for genes which have log-fold-changes different from zero, it tests whether the log2-fold-change is greater than \code{lfc} in absolute value.
\code{glmTreat} is analogous to the TREAT approach developed by McCarthy and Smyth (2009) for microarrays.

Note that the \code{lfc} testing threshold used to define the null hypothesis is not the same as a log2-fold-change cutoff, as the observed log2-fold-change needs to substantially larger than \code{lfc} for the gene to be called as significant.
In practice, modest values for \code{lfc} such as \code{log2(1.1)}, \code{log2(1.2)} or \code{log2(1.5)} are usually the most useful.
In practice, setting \code{lfc=log2(1.2)} or \code{lfc=log2(1.5)} will usually cause most differentially expressed genes to have estimated fold-changes of 2-fold or greater, depending on the sample size and precision of the experiment.

Note also that \code{glmTreat} constructs test statistics using the unshrunk log2-fold-changes(\code{unshrunk.logFC}) rather than the log2-fold-changes that are usually reported (\code{logFC}).
If no shrinkage has been applied to the log-fold-changes, i.e., the glms were fitted with \code{prior.count=0}, then \code{unshrunk.logFC} and \code{logFC} are the same and the former is omitted from the output object.

\code{glmTreat} detects whether \code{glmfit} was produced by \code{glmFit} or \code{glmQLFit}.
In the former case, it conducts a modified likelihood ratio test (LRT) against the fold-change threshold.
In the latter case, it conducts a quasi-likelihood (QL) F-test against the threshold.

If \code{lfc=0}, then \code{glmTreat} is equivalent to \code{glmLRT} or \code{glmQLFTest}, depending on whether likelihood or quasi-likelihood is being used.

\code{glmTreat} with positive \code{lfc} gives larger p-values than would be obtained with \code{lfc=0}.
If \code{null="worst.case"}, then \code{glmTreat} conducts a test closely analogous to the \code{treat} function in the limma package.
This conducts a test if which the null hypothesis puts the true logFC on the boundary of the \code{[-lfc,lfc]} interval closest to the observed logFC.
If \code{null="interval"}, then the null hypotheses assumes an interval of possible values for the true logFC.
This approach is somewhat less conservative.

Note that, unlike other edgeR functions such as \code{glmLRT} and \code{glmQLFTest}, \code{glmTreat} can only accept a single contrast.
If \code{contrast} is a matrix with multiple columns, then only the first column will be used.
}

\note{
\code{glmTreat} was previously called \code{treatDGE} in edgeR versions 3.9.10 and earlier. 
}

\author{Yunshun Chen and Gordon Smyth}

\seealso{
\code{\link{topTags}} displays results from \code{glmTreat}.

\code{\link{treat}} is the corresponding function in the limma package, designed for use with normally distributed log-expression data rather than for negative binomial counts.
}

\references{
McCarthy, D. J., and Smyth, G. K. (2009).
Testing significance relative to a fold-change threshold is a TREAT.
\emph{Bioinformatics} 25, 765-771.
\url{http://bioinformatics.oxfordjournals.org/content/25/6/765}
}

\examples{
ngenes <- 100
n1 <- 3
n2 <- 3
nlibs <- n1+n2
mu <- 100
phi <- 0.1
group <- c(rep(1,n1), rep(2,n2))
design <- model.matrix(~as.factor(group))

### 4-fold change for the first 5 genes
i <- 1:5
fc <- 4
mu <- matrix(mu, ngenes, nlibs)
mu[i, 1:n1] <- mu[i, 1:n1]*fc

counts <- matrix(rnbinom(ngenes*nlibs, mu=mu, size=1/phi), ngenes, nlibs)
d <- DGEList(counts=counts,lib.size=rep(1e6, nlibs), group=group)

gfit <- glmFit(d, design, dispersion=phi)
tr <- glmTreat(gfit, coef=2, lfc=1)
topTags(tr)
}
