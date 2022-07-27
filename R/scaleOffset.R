scaleOffset <- function(y, ...)
UseMethod("scaleOffset")

scaleOffset.default <- function(y, offset, ...)
#	Ensures scale of offsets are consistent with library sizes.
#	Aaron Lun and Yunshun Chen.
#	Created 8 Dec 2016, last modified 1 Jun 2020
{
	if(is.matrix(y)) lib.size <- colSums(y)
	else lib.size <- y

	if(is.matrix(offset)) {
		if (ncol(offset)!=length(lib.size)) {
			stop("'ncol(offset)' should be equal to number of libraries")
		}
		if(inherits(offset,"CompressedMatrix")) {
			if(attr(offset, "repeat.row")) {
				adj <- mean(offset)
			} else if (attr(offset, "repeat.col")) {
				offset <- makeCompressedMatrix(0, dim(offset))
				adj <- 0
			} else {
				adj <- rowMeans(as.matrix(offset))
			}
		} else {
			adj <- rowMeans(offset)
		}
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
#	Created 8 Dec 2016, last modified 29 May 2020
{
	lib.size <- y$samples$lib.size * y$samples$norm.factors
	y$offset <- scaleOffset(lib.size, offset)
	y
}

