\name{read10X}
\alias{read10X}

\title{Read 10X Genomics Files}

\description{Reads 10X Genomics files containing single-cell RNA-seq UMI counts in Matrix Market format.}

\usage{
read10X(mtx = NULL, genes = NULL, barcodes = NULL, path = ".", DGEList = TRUE)
}

\arguments{
\item{mtx}{name of \code{mtx} file containing counts in Matrix Exchange Format. Defaults to \code{matrix.mtx} or \code{matrix.mtx.gz}.}
\item{genes}{name of file containing gene IDs and names. Defaults to \code{features.tsv} or \code{genes.tsv} or gzipped versions of the same.}
\item{barcodes}{optional name of file containing barcodes. Defaults to \code{"barcodes.tsv"} or \code{barcodes.tsv.gz}.}
\item{path}{character string giving the directory containing the files.
Defaults to the current working directory.}
\item{DGEList}{logical. If \code{TRUE}, a \code{DGEList} will be returned, otherwise an unclassed list is returned.}
}

\details{
This function reads output files created by the 10X Genomics Cellranger pipeline, see
\url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices}.
The UMI counts are assembled into an integer matrix in R with accompanying gene IDs and gene symbols.
The results are returned as either a \code{DGEList} or an ordinary list.

The files \code{mtx}, \code{genes} and \code{barcodes} can be provided in either gzipped or unzipped versions.

This function creates an ordinary matrix of counts.
To read the counts instead into a sparse matrix format, the \code{read10xResults} function in the scater package is an alternative.
}

\value{
Either a \code{\link[edgeR:DGEList-class]{DGEList}} object (if \code{DGEList=TRUE}) or an ordinary list with the following components:
\item{counts}{matrix of counts.}
\item{genes}{data.frame counting gene symbols.}
\item{samples}{data.frame containing information about each cell. This will be omitted if \code{barcodes=NULL} and \code{DGEList=FALSE}.}
The only difference between the \code{DGEList} or list formats is that the \code{DGEList} adds some extra columns to the \code{samples} data.frame.
}

\author{Gordon Smyth}

\examples{
\dontrun{
GEO <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2510nnn/GSM2510617/suppl/"
GEOmtx <- paste0(GEO,"GSM2510617_P7-matrix.mtx.gz")
GEOgenes <- paste0(GEO,"GSM2510617_P7-genes.tsv.gz")
download.file(GEOmtx,"matrix.mtx.gz")
download.file(GEOgenes,"genes.tsv.gz")
y <- read10X("matrix.mtx.gz", "genes.tsv.gz")
}}

\seealso{
\code{read10xResults} in the scater package.
}

\concept{Reading data files}
