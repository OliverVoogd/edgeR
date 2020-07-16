effectiveLibSizes <- function(y, log=FALSE, ...)
UseMethod("effectiveLibSizes")

effectiveLibSizes.DGEList <- function(y, log=FALSE, ...)
#	Effective (normalized) library size
#	Gordon Smyth.
#	Created 19 Apr 2020. Last modified.
{
	if(is.null(y$offset)) {
		els <- y$samples$lib.size*y$samples$norm.factors
		if(log) els <- log(els)
	} else {
		els <- y$offset[1,]
		if(!log) els <- exp(els)
	}
	els
}

effectiveLibSizes.DGELRT <- effectiveLibSizes.DGEGLM <- function(y, log=FALSE, ...)
#	Effective (normalized) library size from DGEGLM fit.
#	Gordon Smyth.
#	Created 19 Apr 2020. Last modified.
{
	els <- y$offset[1,]
	if(!log) els <- exp(els)
	els
}

effectiveLibSizes.default <- function(y, log=FALSE, ...)
#	Effective library sizes for a matrix, i.e., just the column sums
#	Gordon Smyth.
#	Created 19 Apr 2020. Last modified.
{
	y <- as.matrix(y)
	els <- colSums(y)
	if(log) els <- log(els)
	els
}
