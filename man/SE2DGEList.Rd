\name{SE2DGEList}
\alias{SE2DGEList}

\title{SummarizedExperiment to DGEList}

\description{Given any \code{SummarizedExperiment} data object, extract basic information needed and convert it into a \code{DGEList} object.}

\usage{
SE2DGEList(object)
}

\arguments{
\item{object}{a \code{SummarizedExperiment} object. Must have \code{counts} in its assay component.}
}

\details{
This function takes a \code{SummarizedExperiment} data object as input.
The \code{counts} of the assay component of the input \code{SummarizedExperiment} data object is extracted and used as the \code{counts} component of the output \code{DGEList} object.
The \code{rowRanges} or \code{rowData} of the input is converted into a \code{data.frame} and used as \code{genes} in the output.
The \code{colData} of the input is also converted into a \code{data.frame} and used as the sample information in the output.
}

\value{
A \code{DGEList} object.
}

\author{Yunshun Chen and Gordon Smyth}

\examples{
\dontrun{
library(SummarizedExperiment)
example(SummarizedExperiment)
y <- SE2DGEList(se)
}}

\concept{edgeR classes}
