cpm <- function(y, ...)
UseMethod("cpm")

cpm.DGEList <- function(y, normalized.lib.sizes=TRUE, log=FALSE, prior.count=2, ...)
#	Counts per million for a DGEList
#	Davis McCarthy and Gordon Smyth.
#	Created 20 June 2011. Last modified 22 October 2020.
{
	lib.size <- y$samples$lib.size
	if(!is.null(y$offset)){
		if( min(y$offset) > max(log(lib.size)) || min(log(lib.size)) > max(y$offset) ) warning("Offset may not reflect library sizes. Scaling offset may be required.")
		lib.size <- NULL
	} else {
		if(normalized.lib.sizes) lib.size <- lib.size*y$samples$norm.factors
	}

	cpm.default(y$counts, lib.size=lib.size, offset=y$offset, log=log, prior.count=prior.count)
}

cpm.SummarizedExperiment <- function(y, normalized.lib.sizes=TRUE, log=FALSE, prior.count=2, ...)
#	Counts per million for a SummarizedExperiment
#	Created 03 April 2020.  Last modified 1 June 2020.
{
	y <- SE2DGEList(y)
	cpm.DGEList(y, normalized.lib.sizes=normalized.lib.sizes, log=log, prior.count=prior.count, ...)
}

cpm.DGELRT <- cpm.DGEGLM <- function(y, log=FALSE, shrunk=TRUE, ...)
#	Fitted counts per million from a fitted model object.
#	Created 19 April 2020.  Last modified 19 April 2020.
{
	if(shrunk) {
		eta <- y$coefficients %*% t(y$design)
	} else {
		eta <- y$unshrunk.coefficients %*% t(y$design)
	}
	(eta + log(1e6)) / log(2)
}

cpm.default <- function(y, lib.size=NULL, offset=NULL, log=FALSE, prior.count=2, ...)
#	Counts per million for a matrix
#	Davis McCarthy and Gordon Smyth.
#	Created 20 June 2011. Last modified 28 May 2020.
{
#	Check y
	y <- as.matrix(y)
	if (any(dim(y)==0L)) {
		return(y)
	}

	if(!is.null(offset)) {
		if(is.matrix(offset)) {
			if(any(dim(offset)!=dim(y))) stop("dimensions are not consistent between counts and offset")
		} else {
			if(length(offset)!=ncol(y)) stop("Length of offset differs from number of libraries")
		}
		if(!is.null(lib.size)) warning("lib.size is ignored in the presence of offset")
		lib.size <- exp(offset)
	} else {
		if(is.null(lib.size)) lib.size <- colSums(y)
	}

	if(!is.double(lib.size)) {
		if(!is.numeric(lib.size)) stop("lib.size must be numeric")
		storage.mode(lib.size) <- "double"
	}

	lib.size <- makeCompressedMatrix(lib.size, dim(y), byrow=TRUE)

	check.range <- suppressWarnings(range(lib.size))
	if (any(is.na(check.range)) || check.range[1] <= 0) {
		stop("library sizes should be finite and non-negative")
	}

#	Calculating in C++ for max efficiency
	if(log) {
		prior.count <- .compressPrior(y, prior.count)
		out <- .Call(.cxx_calculate_cpm_log, y, lib.size, prior.count)
	} else {
		out <- .Call(.cxx_calculate_cpm_raw, y, lib.size)
	}

#	Cleaning up
	dimnames(out) <- dimnames(y)
	out
}
