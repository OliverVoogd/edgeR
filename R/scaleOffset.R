scaleOffset <- function(y, ...)
UseMethod("scaleOffset")

scaleOffset.default <- function(y, offset, ...)
#	Ensures scale of offsets are consistent with library sizes.
#	Aaron Lun and Yunshun Chen.
#	Created 8 Dec 2016.
#   last modified 29 June 2017    
{
	if(is.matrix(y)) lib.size <- colSums(y)
	else lib.size <- y

	if(is.matrix(offset)) {
        if (ncol(offset)!=length(lib.size)) {
            stop("'ncol(offset)' should be equal to number of libraries")
        }
        adj <- rowMeans(offset)       
    } else {
        if (length(offset)!=length(lib.size)) {
            stop("length of 'offset' should be equal to number of libraries")
        }
        adj <- mean(offset)
    }

	mean(log(lib.size)) + offset - adj
}

scaleOffset.DGEList <- function(y, offset, ...) 
#	Ensures scale of offsets are consistent with library sizes.
#	Aaron Lun and Yunshun Chen.
#	Created 8 Dec 2016.
{
    y$offset <- scaleOffset(y$samples$lib.size, offset)
    y
}

