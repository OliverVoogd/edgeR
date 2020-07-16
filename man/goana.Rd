\name{goana.DGELRT}
\alias{goana.DGEExact}
\alias{goana.DGELRT}
\alias{kegga.DGEExact}
\alias{kegga.DGELRT}
\title{Gene Ontology or KEGG Analysis of Differentially Expressed Genes}
\description{
Test for over-representation of gene ontology (GO) terms or KEGG pathways in the up and down differentially expressed genes from a linear model fit.
}
\usage{
\method{goana}{DGELRT}(de, geneid = rownames(de), FDR = 0.05, trend = FALSE, \dots)
\method{kegga}{DGELRT}(de, geneid = rownames(de), FDR = 0.05, trend = FALSE, \dots)
}
\arguments{
  \item{de}{an \code{DGELRT} or \code{DGEExact} object.}
  \item{geneid}{gene IDs. Either a character vector of length \code{nrow(de)} or the name of the column of \code{de$genes} containing the Gene IDs.}
  \item{FDR}{false discovery rate cutoff for differentially expressed genes. Numeric value between 0 and 1.}
  \item{trend}{adjust analysis for gene length or abundance?
  Can be logical, or a numeric vector of covariate values, or the name of the column of \code{de$genes} containing the covariate values.
  If \code{TRUE}, then \code{de$AveLogCPM} is used as the covariate.}
  \item{\dots}{any other arguments are passed to \code{\link{goana.default}} or \code{\link{kegga.default}}.}
}
\details{
\code{goana} performs Gene Ontology enrichment analyses for the up and down differentially expressed genes from a linear model analysis.
\code{kegga} performs the corresponding analysis for KEGG pathways.

The argument \code{de} should be a fitted model object created by \code{glmLRT}, \code{glmTreat}, \code{glmQLFTest} or \code{exactTest}.

For \code{goana}, the gene IDs must be Entrez Gene IDs.
These can be supplied either as row.names of \code{de} or as a column of \code{de$genes}.
In the latter case, the column name containing the Entrez IDs is given by \code{geneid}.
Alternatively, if the Entrez IDs are not part of the \code{de} object, then they can be supplied as a vector argument to \code{geneid}.

For \code{kegga}, gene IDs other than Entrez Gene IDs are supported for some species.
See \code{\link{kegga.default}} for more information.

If \code{trend=FALSE}, the function computes one-sided hypergeometric tests equivalent to Fisher's exact test.

If \code{trend=TRUE} or a covariate is supplied, then a trend is fitted to the differential expression results and the method of Young et al (2010) is used to adjust for this trend.
The adjusted test uses Wallenius' noncentral hypergeometric distribution.
}
\value{
\code{goana} produces a data.frame with a row for each GO term and the following columns:
  \item{Term}{GO term.}
  \item{Ont}{ontology that the GO term belongs to.  Possible values are \code{"BP"}, \code{"CC"} and \code{"MF"}.}
  \item{N}{Number of genes in the GO term.}
  \item{Up}{number of up-regulated differentially expressed genes.}
  \item{Down}{number of down-regulated differentially expressed genes.}
  \item{P.Up}{p-value for over-representation of GO term in up-regulated genes.}
  \item{P.Down}{p-value for over-representation of GO term in down-regulated genes.}
The row names of the data frame give the GO term IDs.

\code{kegga} produces a data.frame as above except that the rownames are KEGG pathway IDs, Term become Path and there is no Ont column.
}

\references{
Chen Y, Lun ATL, and Smyth, GK (2016). From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline.
\emph{F1000Research} 5, 1438. 
\url{https://f1000research.com/articles/5-1438}

  Young, M. D., Wakefield, M. J., Smyth, G. K., Oshlack, A. (2010).
  Gene ontology analysis for RNA-seq: accounting for selection bias.
  \emph{Genome Biology} 11, R14.
  \url{http://genomebiology.com/2010/11/2/R14}
}

\seealso{
\code{\link{goana}}, \code{\link{topGO}}, \code{\link{kegga}}, \code{\link{topKEGG}}
}

\author{Yunshun Chen and Gordon Smyth}

\examples{
\dontrun{
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
go <- goana(lrt, species="Hs")
topGO(go, ont="BP", sort="up")
topGO(go, ont="BP", sort="down")}
}

\concept{Gene set testing}
