\name{readDGE}
\alias{readDGE}

\title{Read and Merge a Set of Files Containing Count Data}

\description{Reads and merges a set of text files containing gene expression counts.}

\usage{readDGE(files, path=NULL, columns=c(1,2), group=NULL, labels=NULL, \dots)}

\arguments{ 
\item{files}{character vector of filenames, or a data.frame of sample information containing a column called \code{files}.}
\item{path}{character string giving the directory containing the files.
Defaults to the current working directory.}
\item{columns}{numeric vector stating which columns of the input files contain the gene names and counts respectively.}
\item{group}{optional vector or factor indicating the experimental group to which each file belongs.}
\item{labels}{character vector giving short names to associate with the files. Defaults to the file names.}
\item{\dots}{other arguments are passed to \code{\link{read.delim}}.}
}

\details{
Each file is assumed to contain digital gene expression data for one genomic sample or count library, with gene identifiers in the first column and counts in the second column.
Gene identifiers are assumed to be unique and not repeated in any one file. 
The function creates a combined table of counts with rows for genes and columns for samples.
A count of zero will be entered for any gene that was not found in any particular sample.

By default, the files are assumed to be tab-delimited and to contain column headings.
Other file formats can be handled by adding arguments to be passed to \code{read.delim}.
For example, use \code{header=FALSE} if there are no column headings and use \code{sep=","} to read a comma-separated file.

Instead of being a vector, the argument \code{files} can be a data.frame containing all the necessary sample information.
In that case, the filenames and group identifiers can be given as columns \code{files} and \code{group} respectively, and the \code{labels} can be given as the row.names of the data.frame.
}

\value{
A \code{\link[edgeR:DGEList-class]{DGEList}} object containing a matrix of counts, with a row for each unique tag found in the input files and a column for each input file.
}

\author{Mark Robinson and Gordon Smyth}

\examples{
#  Read all .txt files from current working directory

\dontrun{files <- dir(pattern="*\\\\.txt$")
RG <- readDGE(files)}
}

\seealso{
See \code{\link{read.delim}} for other possible arguments that can be accepted.

\code{\link{DGEList-class}}, \code{\link{DGEList}}.
}

\concept{Reading data files}
