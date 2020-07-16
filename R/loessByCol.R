loessByCol <- function(y, x=NULL, span=0.5)
# Calls a C++ function to do the dirty work of fitting a degree-0,
# non-robustified loess curve through each column of a matrix.

# C++ version by Aaron Lun, 26 June 2012.  Last modified 6 July 2012.
# Replaces:
# Rcode version by Davis McCarthy, May 2010.
# simpleLoess version by Yunshun Chen, 08 May 2012.
{
	y <- as.matrix(y)
	ntags <- nrow(y)
	if(is.null(x)) x <- 1:ntags

	# Sort by x-values.
	x.order <- order(x)
	y <- y[x.order,,drop=FALSE]
	x <- x[x.order]

	nspan <- min(floor(span*ntags), ntags)
	if(nspan<=1) {
	   fitted <- list(fitted.values=y,leverages=rep(1,ntags))
	   names(fitted$leverages) <- rownames(y)
	   return(fitted)
	}

	# Passing to the compiled code. Note type checking, otherwise the code will complain.
	if (!is.double(y)) storage.mode(y) <- "double"
	if (!is.double(x)) x <- as.double(x)
	fitted <- .Call(.cxx_loess_by_col, x, y, ncol(y), nspan)
   
	# Recover the original order.	
	fitted[[1]][x.order,] <- fitted[[1]]
	fitted[[2]][x.order] <- fitted[[2]]

	# Beautifying.
	names(fitted) <- c("fitted.values", "leverages")
	dimnames(fitted$fitted.values) <- dimnames(y)
	names(fitted$leverages) <- rownames(y)

	fitted
}
