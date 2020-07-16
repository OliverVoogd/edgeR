modelMatrixMeth <- function(object, ...)
#	Create expanded model matrix (aka design matrix) for edgeR analysis of BS-seq methylation counts.
#	Gordon Smyth
#	Created 26 Dec 2017.
{
#	If object isn't already a design matrix, create the design matrix from treatments at the sample level
	if(is.matrix(object)) {
		design.treatments <- object
	} else {
		design.treatments <- model.matrix(object, ...)
	}

#	Design matrix for samples
#	Allow for two observations (methylated and unmethylated counts) for each sample
	nsamples <- nrow(design.treatments)
	Sample <- gl(nsamples,2L)
	design.samples <- model.matrix(~0+Sample)

#	Expand design.treatments for meth & unmeth counts
	design.treatments <- design.treatments[as.integer(Sample),]

#	Combine with treatments applying to methylation effects
	Methylation <- gl(2L,1L,2L*nsamples)
	cbind(design.samples, (Methylation==1L) * design.treatments)
}
