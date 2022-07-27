adjustedProfileLik <- function(dispersion, y, design, offset, weights=NULL, adjust=TRUE, start=NULL, get.coef=FALSE)
#	Tagwise Cox-Reid adjusted profile log-likelihoods for the dispersion.
#	dispersion can be a scalar or a tagwise vector.
#	Computationally, dispersion can also be a matrix, but the apl is still computed tagwise.
#	y is a matrix: rows are genes/tags/transcripts, columns are samples/libraries.
#	offset is a matrix of the same dimensions as y.

#	The weights argument was added by Xiaobei Zhou 20 March 2013,
#	but the log NB probabilities were incorrectly multiplied by the weights.
#	This is fixed 1 March 2018 with a more rigorous interpretation of weights
#	in terms of averages.

#	Yunshun Chen, Gordon Smyth, Aaron Lun
#	Created June 2010. Last modified 22 May 2020.
{
#	Checking counts
	if (!is.numeric(y)) stop("counts must be numeric")
	y <- as.matrix(y)

#	Checking offsets
	offset <- .compressOffsets(y, offset=offset)

#	Checking dispersion
	dispersion <- .compressDispersions(y, dispersion)

#	Checking weights
	weights <- .compressWeights(y, weights)
	  
#	Fit tagwise linear models
	fit <- glmFit(y,design=design,dispersion=dispersion,offset=offset,prior.count=0,weights=weights,start=start)
	mu <- fit$fitted.values

#	Check other inputs to C++ code
	adjust <- as.logical(adjust)
	if (!is.double(design)) storage.mode(design) <- "double"

#	Compute adjusted log-likelihood
	apl <- .Call(.cxx_compute_apl, y, mu, dispersion, weights, adjust, design)

#	Deciding what to return.
	if (get.coef) { 
		return(list(apl=apl, beta=fit$coefficients))
	} else {
		return(apl)
	}
}

