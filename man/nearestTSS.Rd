\name{nearestTSS}
\alias{nearestTSS}
\title{Find Nearest Transcriptional Start Site}

\description{
Find nearest TSS and distance to nearest TSS for a vector of chromosome loci.
}

\usage{
nearestTSS(chr, locus, species="Hs")
}

\arguments{
  \item{chr}{character vector of chromosome names.}
  \item{locus}{integer or numeric vector of genomic loci, of same length as \code{chr}.}
  \item{species}{character string specifying the species.
  Possible values include \code{"Hs"} (human), \code{"Mm"} (mouse), \code{"Rn"} (rat), \code{"Dm"} (fly) or \code{"Pt"} (chimpanzee), but other values are possible if the corresponding organism package is available.
  See \code{\link{alias2Symbol}} for other possible values.}
}

\details{
This function takes a series of genomic loci, defined by a vector of chromosome names and a vector of genomic positions within the chromosomes,
and finds the nearest transcriptional start site (TSS) for each locus.
The chromosome names can be in the format \code{"1","2","X"} or can be \code{"chr1","chr2","chrX"}.

This function uses the Bioconductor organism package named "org.XX.eg.db" where XX is \code{species}.
}

\value{
A data.frame with components:
  \item{gene_id}{character vector giving the Entrez Gene ID of the nearest TSS for each element of \code{chr} and \code{locus}.}
  \item{symbol}{character vector of gene symbols.}
  \item{strand}{character vector with \code{"+"} for positive strand genes and \code{"-"} for negative strand genes.}
  \item{tss}{integer vector giving TSS.}
  \item{width}{integer vector giving genomic width of the gene.}
  \item{distance}{integer vector giving distance to nearest TSS. Positive values means that the TSS is downstream of the locus, negative values means that it is upstream. Gene body loci will therefore have negative distances and promotor loci will have positive.}
}

\author{Gordon Smyth}

\seealso{
\code{\link{nearestReftoX}}
}

\examples{
nearestTSS(chr = c("1","1"), locus = c(1000000,2000000))
}
