addPriorCount <- function(y, lib.size=NULL, offset=NULL, prior.count=1) 
# Add library size-adjusted prior counts to values of 'y'.
# Also add twice the adjusted prior to the library sizes, 
# which are provided as log-transformed values in 'offset'.
#
# written by Aaron Lun
# created 26 September 2016
# last modified 4 Nov 2018    
{
#	Check y
	y <- as.matrix(y)
	if (!is.numeric(y)) stop('count matrix must be numeric')

#	Check prior.count
	prior.count <- .compressPrior(y, prior.count)

#	Check lib.size and offset.
#	If offsets are provided, they must have a similar average to log(lib.size)
#	for the results to be meaningful as logCPM values
	offset <- .compressOffsets(y, lib.size=lib.size, offset=offset)

#	Adding the prior count.
	out <- .Call(.cxx_add_prior_count, y, offset, prior.count)
	names(out) <- c("y", "offset")
	dimnames(out$y) <- dimnames(y)
	out$offset <- makeCompressedMatrix(out$offset, dim(y), byrow=TRUE)
	return(out)
}

