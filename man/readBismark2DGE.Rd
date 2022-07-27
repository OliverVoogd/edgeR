\name{readBismark2DGE}
\alias{readBismark2DGE}

\title{Read Bismark Coverage Files}

\description{Read Bismark coverage files containing methylated and unmethylated read counts for CpG loci and create DGEList.}

\usage{readBismark2DGE(files, sample.names=NULL, readr=TRUE, verbose=TRUE)}

\arguments{ 
\item{files}{character vector of file names.}
\item{sample.names}{character vector of sample names. If \code{NULL}, sample names will be extracted from the file names.}
\item{readr}{logical. If \code{TRUE}, \code{readr} package is used to read the coverage files, otherwise \code{read.delim} is used.}
\item{verbose}{logical. If \code{TRUE}, read progress messages are send to standard output.}
}

\details{
This function reads tab-delimited coverage files output by Bismark software.
Counts from multiple files are collated into a \code{\link[edgeR:DGEList-class]{DGEList}} object.
}

\note{
This function represents genomic loci as integers, so the largest locus position must be less than the maximum integer in R (about \code{2e9}).
The number of chromosomes times the largest locus position must be less than \code{1e16}.
}

\value{
A \code{\link[edgeR:DGEList-class]{DGEList}} object with a row for each unique genomic loci found in the files and two columns (containing methylated and unmethylated counts) for each sample.
}

\author{Gordon Smyth}

\references{
Chen, Y, Pal, B, Visvader, JE, Smyth, GK (2017).
Differential methylation analysis of reduced representation bisulfite sequencing experiments using edgeR.
\emph{F1000Research} 6, 2055.
\url{https://f1000research.com/articles/6-2055} 
}

\concept{Reading data files}
\concept{Differential methylation}
