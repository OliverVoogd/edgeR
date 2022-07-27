\name{topTags}
\alias{topTags}
\alias{TopTags-class}
\alias{show,TopTags-method}

\title{Table of the Top Differentially Expressed Genes/Tags}

\description{Extracts the most differentially expressed genes (or sequence tags) from a test object, ranked either by p-value or by absolute log-fold-change.}

\usage{
topTags(object, n = 10, adjust.method = "BH", sort.by = "PValue", p.value = 1)
}

\arguments{ 
\item{object}{a \code{\link[edgeR:DGEList-class]{DGEExact}} or \code{\link[edgeR:DGELRT-class]{DGELRT}} object containing test statistics and p-values.
Usually created by \code{exactTest}, \code{glmLRT}, \code{glmTreat} or \code{glmQLFTest}.}

\item{n}{integer, maximum number of genes/tags to return.}

\item{adjust.method}{character string specifying the method used to adjust p-values for multiple testing. See \code{\link{p.adjust}} for possible values.}

\item{sort.by}{character string specifying the sort method. Possibilities are \code{"PValue"} for p-value, \code{"logFC"} for absolute log-fold change or \code{"none"} for no sorting.}

\item{p.value}{numeric cutoff value for adjusted p-values. Only tags with adjusted p-values equal or lower than specified are returned.}
}

\details{
This function is closely analogous to the \code{\link{topTable}} function in the limma package.
It accepts a test statistic object created by any of the edgeR functions \code{exactTest}, \code{glmLRT}, \code{glmTreat} or \code{glmQLFTest} and extracts a readable data.frame of the most differentially expressed genes.
The data.frame collates the annotation and differential expression statistics for the top genes.
The data.frame is wrapped in a \code{TopTags} output object that records the test statistic used and the multiple testing  adjustment method.

\code{TopTags} objects will return dimensions and hence functions such as \code{dim}, \code{nrow} or \code{ncol} are defined on them.
\code{TopTags} objects also have a \code{show} method so that printing produces a compact summary of their contents.

\code{topTags} permits ranking by fold-change but the authors do not recommend fold-change ranking or fold-change cutoffs for routine RNA-seq analysis.
The p-value ranking is intended to more biologically meaningful, especially if the p-values were computed using \code{glmTreat}.
}

\value{
An object of class \code{TopTags}, which is a list-based class with the following components:
\item{table}{a data.frame containing differential expression results for the top genes in sorted order.
The number of rows is the smaller of \code{n} and the number of genes with adjusted p-value less than or equal to \code{p.value}.
The data.frame includes all the annotation columns from \code{object$genes} and all statistic columns from \code{object$table} plus one of:
\tabular{rl}{
\code{FDR}: \tab false discovery rate (only when \code{adjust.method} is \code{"BH"}, \code{"BY"} or \code{"fdr"})\cr
\code{FWER}: \tab family-wise error rate (only when \code{adjust.method} is \code{"holm"}, \code{"hochberg"}, \code{"hommel"} or \code{"bonferroni"}).
}
}
\item{adjust.method}{character string specifying the method used to adjust p-values for multiple testing, same as input argument.}
\item{comparison}{character vector giving the names of the two groups being compared (for \code{DGEExact} objects) or the glm contrast being tested (for \code{DGELRT} objects).}
\item{test}{character string stating the name of the test.}
}

\note{
The terms `tag' and `gene' are used synonymously on this page and refer to the rows of \code{object}.
In general, the rows might be genes, sequence tags, transcripts, exons or whatever type of genomic feature is appropriate for the analysis at hand.
}

\author{Mark Robinson, Davis McCarthy, Yunshun Chen, Gordon Smyth}

\references{
Chen Y, Lun ATL, and Smyth, GK (2016).
From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline.
\emph{F1000Research} 5, 1438.
\url{http://f1000research.com/articles/5-1438}

McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
\emph{Nucleic Acids Research} 40, 4288-4297.
\doi{10.1093/nar/gks042}

Robinson MD, Smyth GK (2008). Small-sample estimation of negative binomial dispersion, with applications to SAGE data. \emph{Biostatistics} 9, 321-332.

Robinson MD, Smyth GK (2007). Moderated statistical tests for assessing differences in tag abundance. \emph{Bioinformatics} 23, 2881-2887.
}

\seealso{
\code{\link{exactTest}}, \code{\link{glmLRT}}, \code{\link{glmTreat}}, \code{\link{glmQLFTest}}, \code{\link{dim.TopTags}}, \code{\link{p.adjust}}.
}

\examples{
# generate raw counts from NB, create list object
y <- matrix(rnbinom(80,size=1,mu=10),nrow=20)
d <- DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
rownames(d$counts) <- paste("gene",1:nrow(d$counts),sep=".")

# estimate common dispersion and find differences in expression
# here we demonstrate the 'exact' methods, but the use of topTags is
# the same for a GLM analysis
d <- estimateCommonDisp(d)
de <- exactTest(d)

# look at top 10
topTags(de)
# Can specify how many genes to view
tp <- topTags(de, n=15)
# Here we view top 15
tp
# Or order by fold change instead
topTags(de,sort.by="logFC")
}
