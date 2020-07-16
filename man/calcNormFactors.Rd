\name{calcNormFactors}
\alias{calcNormFactors}
\alias{calcNormFactors.DGEList}
\alias{calcNormFactors.SummarizedExperiment}
\alias{calcNormFactors.default}

\title{Calculate Normalization Factors to Align Columns of a Count Matrix}
\description{
Calculate normalization factors to scale the raw library sizes.
}
\usage{
\method{calcNormFactors}{DGEList}(object,
                method=c("TMM","TMMwsp","RLE","upperquartile","none"),
                refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE,
                Acutoff=-1e10, p=0.75, \dots)
\method{calcNormFactors}{SummarizedExperiment}(object,
                method=c("TMM","TMMwsp","RLE","upperquartile","none"),
                refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE,
                Acutoff=-1e10, p=0.75, \dots)
\method{calcNormFactors}{default}(object, lib.size=NULL,
                method=c("TMM","TMMwsp","RLE","upperquartile","none"),
                refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE,
                Acutoff=-1e10, p=0.75, \dots)
}

\arguments{
  \item{object}{a \code{matrix} of raw (read) counts, or a \code{DGEList} object, or a \code{SummarizedExperiment} object}
  \item{lib.size}{numeric vector of library sizes of the \code{object}.}
  \item{method}{normalization method to be used}
  \item{refColumn}{column to use as reference for \code{method="TMM"}. Can be a column number or a numeric vector of length \code{nrow(object)}.}
  \item{logratioTrim}{amount of trim to use on log-ratios ("M" values) for \code{method="TMM"}}
  \item{sumTrim}{amount of trim to use on the combined absolute levels ("A" values) for \code{method="TMM"}}
  \item{doWeighting}{logical, whether to compute (asymptotic binomial precision) weights for \code{method="TMM"}}
  \item{Acutoff}{cutoff on "A" values to use before trimming for \code{method="TMM"}}
  \item{p}{percentile (between 0 and 1) of the counts that is aligned when \code{method="upperquartile"}}
  \item{\dots}{other arguments that are not currently used.}
}

\details{
This function computes scaling factors to convert observed library sizes into effective library sizes.
The effective library sizes for use in downstream analysis are \code{lib.size * norm.factors} where \code{lib.size} contains the original library sizes and \code{norm.factors} is the output from this function.

The TMM method implements the trimmed mean of M-values proposed by Robinson and Oshlack (2010).
By default, the M-values are weighted according to inverse variances, as computed by the delta method for logarithms of binomial random variables.
If \code{refColumn} is unspecified, then the library whose upper quartile is closest to the mean upper quartile is used.

\code{method="TMMwsp"} stands for "TMM with singleton pairing".
This is a variant of TMM that is intended to perform better for data with a high proportion of zeros.
In the TMM method, genes that have zero count in either library are ignored when comparing pairs of libraries.
In the TMMwsp method, the positive counts from such genes are reused to increase the number of features by which the libraries are compared.
The singleton positive counts are paired up between the libraries in decreasing order of size and then a slightly modified TMM method is applied to the re-ordered libraries.

RLE is the scaling factor method proposed by Anders and Huber (2010).
We call it "relative log expression", as median library is calculated from the geometric mean of all columns and the median ratio of each sample to the median library is taken as the scale factor.

The upperquartile method is the upper-quartile normalization method of Bullard et al (2010), in which the scale factors are calculated from the 75\% quantile of the counts for each library, after removing genes that are zero in all libraries.
This idea is generalized here to allow scaling by any quantile of the distributions.

If \code{method="none"}, then the normalization factors are set to 1.

For symmetry, normalization factors are adjusted to multiply to 1.
Rows that have zero counts for all columns are removed before normalization factors are computed, so that rows with all zero counts do not affect the estimated factors.
}

\value{
If \code{object} is a \code{matrix}, the output is a vector with length \code{ncol(object)} giving the relative normalization factors.
If \code{object} is a \code{DGEList}, then it is returned as output with the relative normalization factors in \code{object$samples$norm.factors}.
If \code{object} is a \code{SummarizedExperiment}, then it is converted in to a \code{DGEList} object with the relative normalization factors in \code{object$samples$norm.factors}.
}

\author{Mark Robinson, Gordon Smyth}

\references{
Anders, S, Huber, W (2010).
Differential expression analysis for sequence count data
\emph{Genome Biology} 11, R106.

Bullard JH, Purdom E, Hansen KD, Dudoit S. (2010)
Evaluation of statistical methods for normalization and differential expression in mRNA-Seq experiments.
\emph{BMC Bioinformatics} 11, 94.

Robinson MD, Oshlack A (2010).
A scaling normalization method for differential expression analysis of RNA-seq data.
\emph{Genome Biology} 11, R25.
}

\examples{
y <- matrix( rpois(1000, lambda=5), nrow=200 )
calcNormFactors(y)

# The TMMwsp and TMM methods ignore genes with largest fold-changes:
y <- cbind(1,c(1,1,1,1,1,1,1,1,1,100))
calcNormFactors(y, lib.size=c(1e6,1e6))

# calcNormFactors makes the fold-changes for the majority of genes as close to 1 as possible:
# In this example, normalizing the library sizes makes most of the CPMs equal in the two samples:
dge <- DGEList(counts=y)
cpm(dge)
dge <- calcNormFactors(dge)
cpm(dge)
}

\concept{Normalization}
