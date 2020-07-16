\name{edgeR-package}
\docType{package}
\alias{edgeR}
\alias{edgeR-package}
\title{Empirical analysis of digital gene expression data in R}
\description{
edgeR is a package for the analysis of digital gene expression data arising from RNA sequencing technologies such as SAGE, CAGE, Tag-seq or RNA-seq, with emphasis on testing for differential expression.
It can also be used for other sequencing technologies from which read counts are produced, such as ChIP-seq, Hi-C or CRISPR.

Particular strengths of the package include the ability to estimate biological variation between replicate libraries, and to conduct exact tests of significance which are suitable for small counts.
The package is able to make use of even minimal numbers of replicates.

The supplied counts are assumed to be those of genes in a RNA-seq experiment.
However, counts can be supplied for any genomic feature of interest, e.g., tags, transcripts, exons, or even arbitrary intervals of the genome.

An extensive User's Guide is available, and can be opened by typing \code{edgeRUsersGuide()} at the R prompt.
Detailed help pages are also provided for each individual function.

The edgeR package implements original statistical methodology described in the publications below.
}

\author{
Yunshun Chen, Aaron TL Lun, Davis J McCarthy, Matthew E Ritchie, Belinda Phipson, Yifang Hu, Xiaobei Zhou, Mark D Robinson, Gordon K Smyth
}

\references{
  Robinson MD and Smyth GK (2007). Moderated statistical tests for
  assessing differences in tag abundance. \emph{Bioinformatics} 23, 2881-2887

  Robinson MD and Smyth GK (2008). Small-sample estimation of negative
  binomial dispersion, with applications to SAGE data. \emph{Biostatistics},
  9, 321-332

  Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor
  package for differential expression analysis of digital gene
  expression data. \emph{Bioinformatics} 26, 139-140

  Robinson, MD, Oshlack, A (2010). A scaling normalization method for differential expression analysis of RNA-seq data. \emph{Genome Biology}, 11(3), R25.

  McCarthy, DJ, Chen, Y, Smyth, GK (2012).
  Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
  \emph{Nucleic Acids Research} 40, 4288-4297.

  Anders, S, McCarthy, DJ, Chen, Y, Okoniewski, M, Smyth, GK, Huber, W, and Robinson, MD (2013). Count-based differential expression analysis of RNA sequencing data using R and Bioconductor.
  \emph{Nature Protocols} 8, 1765-1786.

  Zhou, X, Lindsay, H, Robinson, MD (2014). Robustly detecting differential expression in RNA sequencing data using observation weights.
  \emph{Nucleic Acids Research} 42(11), e91.

  Chen, Y, Lun, ATL, and Smyth, GK (2014).
  Differential expression analysis of complex RNA-seq experiments using edgeR.
  In: \emph{Statistical Analysis of Next Generation Sequence Data},
  Somnath Datta and Daniel S Nettleton (eds), Springer, New York, pages 51-74.
  \url{http://www.statsci.org/smyth/pubs/edgeRChapterPreprint.pdf}

  Dai Z, Sheridan, JM, Gearing, LJ, Moore, DL, Su, S, Wormald, S, Wilcox, S, O'Connor, L, Dickins, RA, Blewitt, ME, Ritchie, ME (2014).
  edgeR: a versatile tool for the analysis of shRNA-seq and CRISPR-Cas9 genetic screens.
  \emph{F1000Research} 3, 95.
  \url{http://f1000research.com/articles/3-95}

  Lun, ATL, Chen, Y, and Smyth, GK (2016).
  It's DE-licious: a recipe for differential expression analyses of RNA-seq experiments using quasi-likelihood methods in edgeR.
  \emph{Methods in Molecular Biology} 1418, 391-416.
  \url{http://www.statsci.org/smyth/pubs/QLedgeRPreprint.pdf"} (Preprint 8 April 2015)

  Chen Y, Lun ATL, and Smyth, GK (2016).
  From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline.
  \emph{F1000Research} 5, 1438.
  \url{http://f1000research.com/articles/5-1438}

  Lun, AT, Smyth, GK (2017). No counts, no variance: allowing for loss of degrees of freedom when assessing biological variability from RNA-seq data.
  \emph{Statistical Applications in Genetics and Molecular Biology} 16(2), 83-93.

  Chen, Y, Pal, B, Visvader, JE, Smyth, GK (2017).
  Differential methylation analysis of reduced representation bisulfite sequencing experiments using edgeR.
  \emph{F1000Research} 6, 2055.
}

\concept{Documentation}
