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
  Possible values are \code{"Hs"} (human hg38), \code{"Mm"} (mouse mm10), \code{"Rn"} (rat), \code{"Dm"} (fly), \code{"Dr"} (zebra fish), \code{"Ce"} (worm), \code{"Bt"} (bovine), \code{"Gg"} (chicken), \code{"Mmu"} (rhesus), \code{"Cf"} (canine) or \code{"Pt"} (chimpanzee).}
}

\details{
This function takes a series of genomic loci, defined by a vector of chromosome names and a vector of genomic positions within the chromosomes,
and finds the nearest transcriptional start site (TSS) for each locus.
The chromosome names can be in the format \code{"1","2","X"} or can be \code{"chr1","chr2","chrX"}.

For genes with more than one annotated TSS, only the most 5' (upstream) of the alternative TSS is reported.

This function uses the Bioconductor organism package named "org.XX.eg.db" where XX is \code{species}.
Note that each organism package supports only a particular build of the genome for that species.
For human (\code{species="Hs"}, the results are for the hg38 genome build.
For mouse (\code{species="Mm"}), the results are for the mm10 genome build.
}

\value{
A data.frame with the following columns:
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
