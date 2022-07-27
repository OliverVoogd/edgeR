\name{spliceVariants}
\alias{spliceVariants}

\title{Identify Genes with Splice Variants}

\description{Identify genes exhibiting evidence for splice variants (alternative exon usage/transcript isoforms) from exon-level count data using negative binomial generalized linear models.}

\usage{
spliceVariants(y, geneID, dispersion=NULL, group=NULL, estimate.genewise.disp=TRUE,
               trace=FALSE)
}
\arguments{ 

\item{y}{either a matrix of exon-level counts or a \code{DGEList} object with (at least) elements \code{counts} (table of counts summarized at the exon level) and \code{samples} (data frame containing information about experimental group, library size and normalization factor for the library size). Each row of \code{y} should represent one exon.}

\item{geneID}{vector of length equal to the number of rows of \code{y}, which provides the gene identifier for each exon in \code{y}. These identifiers are used to group the relevant exons into genes for the gene-level analysis of splice variation.}

\item{dispersion}{scalar (in future a vector will also be allowed) supplying the negative binomial dispersion parameter to be used in the negative binomial generalized linear model.}

\item{group}{factor supplying the experimental group/condition to which each sample (column of \code{y}) belongs. If \code{NULL} (default) the function will try to extract if from \code{y}, which only works if \code{y} is a \code{DGEList} object.}

\item{estimate.genewise.disp}{logical, should genewise dispersions (as opposed to a common dispersion value) be computed if the \code{dispersion} argument is \code{NULL}?}

\item{trace}{logical, whether or not verbose comments should be printed as function is run. Default is \code{FALSE}.}
}

\value{\code{spliceVariants} returns a \code{DGEExact} object, which contains a table of results for the test of differential splicing between experimental groups (alternative exon usage), a data frame containing the gene identifiers for which results were obtained and the dispersion estimate(s) used in the statistical models and testing.}

\details{
This function can be used to identify genes showing evidence  of splice variation (i.e. alternative splicing, alternative exon usage, transcript isoforms). A negative binomial generalized linear model is used to assess evidence, for each gene, given the counts for the exons for each gene, by fitting a model with an interaction between exon and experimental group and comparing this model (using a likelihood ratio test) to a null model which does not contain the interaction. Genes that show significant evidence for an interaction between exon and experimental group by definition show evidence for splice variation, as this indicates that the observed differences between the exon counts between the different experimental groups cannot be explained by consistent differential expression of the gene across all exons. The function \code{topTags} can be used to display the results of \code{spliceVariants} with genes ranked by evidence for splice variation.
}

\author{Davis McCarthy, Gordon Smyth}

\examples{
# generate exon counts from NB, create list object
y<-matrix(rnbinom(40,size=1,mu=10),nrow=10)
d<-DGEList(counts=y,group=rep(1:2,each=2))
genes <- rep(c("gene.1","gene.2"), each=5)
disp <- 0.2
spliceVariants(d, genes, disp)
}

\seealso{
\code{\link{estimateExonGenewiseDisp}} for more information about estimating genewise dispersion values from exon-level counts. \code{\link{DGEList}} for more information about the \code{DGEList} class. \code{\link{topTags}} for more information on displaying ranked results from \code{spliceVariants}. \code{\link{estimateCommonDisp}} and related functions for estimating the dispersion parameter for the negative binomial model.
}

\concept{Differential exon usage}
