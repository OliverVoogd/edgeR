\name{normalizeChIPtoInput}
\alias{normalizeChIPtoInput}
\alias{calcNormOffsetsforChIP}

\title{Normalize ChIP-Seq Read Counts to Input and Test for Enrichment}

\description{
Normalize ChIP-Seq read counts to input control values, then test for significant enrichment relative to the control. 
}

\usage{
normalizeChIPtoInput(input, response, dispersion=0.01, niter=6, loss="p", plot=FALSE,
                     verbose=FALSE, \dots)
calcNormOffsetsforChIP(input, response, dispersion=0.01, niter=6, loss="p", plot=FALSE,
                       verbose=FALSE, \dots)
}

\arguments{
  \item{input}{numeric vector of non-negative input values, not necessarily integer.}
  \item{response}{vector of non-negative integer counts of some ChIP-Seq mark for each gene or other genomic feature.}
  \item{dispersion}{negative binomial dispersion, must be positive.}
  \item{niter}{number of iterations.}
  \item{loss}{loss function to be used when fitting the response counts to the input: \code{"p"} for cumulative probabilities or \code{"z"} for z-value.}
  \item{plot}{if \code{TRUE}, a plot of the fit is produced.}
  \item{verbose}{if \code{TRUE}, working estimates from each iteration are output.}
  \item{\ldots}{other arguments are passed to the \code{plot} function.}
}

\details{
\code{normalizeChIPtoInput} identifies significant enrichment for a ChIP-Seq mark relative to input values.
The ChIP-Seq mark might be for example transcriptional factor binding or an epigenetic mark.
The function works on the data from one sample.
Replicate libraries are not explicitly accounted for; this function can either be run on each sample individually or on a pooled of replicates.

ChIP-Seq counts are assumed to be summarized by gene or similar genomic feature of interest.

This function makes the assumption that a non-negligible proportion of the genes, say 25\% or more, are not truly marked by the ChIP-Seq feature of interest.
Unmarked genes are further assumed to have counts at a background level proportional to the input.
The function aligns the counts to the input so that the counts for the unmarked genes behave like a random sample.
The function estimates the proportion of marked genes, and removes marked genes from the fitting process.
For this purpose, marked genes are those with a Holm-adjusted mid-p-value less than 0.5.

When \code{plot=TRUE}, the genes shown in red are the marked genes (with Holm mid-p-value < 0.5) that have been removed as probably enriched during the fitting process.
The normalization line has been fitted to the non-marked genes plotted in black.

The read counts are treated as negative binomial.
The dispersion parameter is not estimated from the data; instead a reasonable value is assumed to be given.

\code{calcNormOffsetsforChIP} returns a numeric matrix of offsets, ready for linear modelling.
}

\value{
\code{normalizeChIPtoInput} returns a list with components
	\item{p.value}{numeric vector of p-values for enrichment.}
	\item{scaling.factor}{factor by which input is scaled to align with response counts for unmarked genes.}
	\item{prop.enriched}{proportion of marked genes, as internally estimated}

\code{calcNormOffsetsforChIP} returns a numeric matrix of offsets.
}

\author{Gordon Smyth}

\concept{Normalization}
