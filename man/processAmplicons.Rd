\name{processAmplicons}
\alias{processAmplicons}

\title{Process raw fastq data from pooled genetic sequencing screens}

\description{
Given a list of sample-specific index (barcode) sequences and hairpin/sgRNA-specific sequences from an amplicon sequencing screen, generate a DGEList of counts from the raw fastq file/(s) containing the sequence reads. The position of the index sequences and hairpin/sgRNA sequences is considered variable, with the hairpin/sgRNA sequences assumed to be located after the index sequences in the read.
}

\usage{
processAmplicons(readfile, readfile2=NULL, barcodefile, hairpinfile,
                    allowMismatch=FALSE, barcodeMismatchBase=1,
                    hairpinMismatchBase=2, dualIndexForwardRead=FALSE,
                    verbose=FALSE, barcodesInHeader=FALSE,
                    plotPositions=FALSE)
}

\arguments{
\item{readfile}{character vector giving one or more fastq filenames}
\item{readfile2}{character vector giving one or more fastq filenames for reverse read, default to NULL}
\item{barcodefile}{filename containing sample-specific barcode ids and sequences}
\item{hairpinfile}{filename containing hairpin/sgRNA-specific ids and sequences}
\item{allowMismatch}{logical, indicates whether sequence mismatch is allowed}
\item{barcodeMismatchBase}{numeric value of maximum number of base sequence mismatches allowed in a barcode sequence when \code{allowMismatch} is \code{TRUE}}
\item{hairpinMismatchBase}{numeric value of maximum number of base sequence mismatches allowed in a hairpin/sgRNA sequence when \code{allowMismatch} is \code{TRUE}}
\item{dualIndexForwardRead}{logical, indicates if forward reads contains a second barcode sequence (must be present in
\code{barcodefile}) which should be matched}
\item{verbose}{if \code{TRUE}, output program progress}
\item{barcodesInHeader}{logical, indicates if barcode sequences should be matched in the header (sequence identifier) of each read (i.e. the first of every group of four lines in the fastq files)}
\item{plotPositions}{logical, indicates if a density plot displaying the position of each barcode and hairpin/sgRNA sequence in the reads should be created. If \code{dualIndexForwardRead} is \code{TRUE} or \code{readfile2} is not \code{NULL}, plotPositions will generate two density plots, side by side, indicating the positions of the first barcodes and hairpins in the first plot, and second barcodes in the second.}
}

\value{Returns a \code{\link[edgeR:DGEList-class]{DGEList}} object with following components:
	\item{counts}{read count matrix tallying up the number of reads with particular barcode and hairpin/sgRNA matches. Each row is a hairpin/sgRNA and each column is a sample}
	\item{genes}{In this case, hairpin/sgRNA-specific information (ID, sequences, corresponding target gene) may be recorded in this data.frame}
	\item{lib.size}{auto-calculated column sum of the counts matrix}
}

\details{
The \code{processAmplicons} function allows for hairpins/sgRNAs/sample index sequences to be in variable positions within each read.

The input barcode file and hairpin/sgRNA files are tab-separated text files with at least two columns (named 'ID' and 'Sequences') containing the sample or hairpin/sgRNA ids and a second column indicating the sample index or hairpin/sgRNA sequences to be matched.
If \code{dualIndexForwardRead} is \code{TRUE}, a third column 'Sequences2' is expected in the barcode file.
If \code{readfile2} is specified, another column 'SequencesReverse' is expected in the barcode file.
The barcode file may also contain a 'group' column that indicates which experimental group a sample belongs to.
Additional columns in each file will be included in the respective \code{$samples} or \code{$genes} data.frames of the final code{\link[edgeR:DGEList-class]{DGEList}} object.
These files, along with the fastq file/(s) are assumed to be in the current working directory.

To compute the count matrix, matching to the given barcodes and hairpins/sgRNAs is conducted in two rounds.
The first round looks for an exact sequence match for the given barcode sequences and hairpin/sgRNA sequences through the entire read, returning the first match found.
If a match isn't found, the program performs a second round of matching which allows for sequence mismatches if \code{allowMismatch} is set to \code{TRUE}.
The maximum number of mismatch bases in barcode and hairpin/sgRNA are specified by the parameters \code{barcodeMismatchBase} and \code{hairpinMismatchBase} respectively.

The program outputs a \code{\link[edgeR:DGEList-class]{DGEList}} object, with a count matrix indicating the number of times each barcode and hairpin/sgRNA combination could be matched in reads from input fastq file(s).

For further examples and data, refer to the case studies available from \url{http://bioinf.wehi.edu.au/shRNAseq}.
}

\note{
This function replaced the earlier function \code{processHairpinReads} in edgeR 3.7.17.

This function replaces the previous \code{processAmplicons} function, which expected the sequences in the fastq files to have a fixed structure (as per Figure 1A of Dai et al., 2014). This function can be used, and is intended for, reads where hairpins/sgRNAs/sample index sequences can be in variable positions within each read. When \code{plotPositions=TRUE} a density plot of the  match positions is created to allow the user to assess whether they occur in the expected postions.
}

\author{Oliver Voogd, Zhiyin Dai, Shian Su and Matthew Ritchie}

\references{
Dai Z, Sheridan JM,  Gearing, LJ, Moore, DL, Su, S, Wormald, S, Wilcox, S, O'Connor, L, Dickins, RA, Blewitt, ME, Ritchie, ME(2014).
edgeR: a versatile tool for the analysis of shRNA-seq and CRISPR-Cas9 genetic screens.
\emph{F1000Research} 3, 95.
\url{http://f1000research.com/articles/3-95}
}
