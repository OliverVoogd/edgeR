mglmOneGroup <- function(y,dispersion=0,offset=0,weights=NULL,coef.start=NULL,maxit=50,tol=1e-10,verbose=FALSE)
#	Fit single-group negative-binomial glm
#	Aaron Lun and Gordon Smyth
#	18 Aug 2010. Last modified 9 July 2017. 
{
#	Check y
	y <- as.matrix(y)
	if(!is.numeric(y)) stop("y is non-numeric")
	.isAllZero(y)

#	Check dispersion
	dispersion <- .compressDispersions(y, dispersion)

#	Check offset
	offset <- .compressOffsets(y, offset=offset)

#	Check starting values
	if (is.null(coef.start)) coef.start <- NA_real_
	if (!is.double(coef.start)) storage.mode(coef.start) <- "double"
	coef.start <- rep(coef.start, length.out=nrow(y))

#	Check weights
	weights <- .compressWeights(y, weights)

#	Fisher scoring iteration.
	output <- .Call(.cxx_fit_one_group, y, offset, dispersion, weights, maxit, tol, coef.start)

#	Convergence achieved for all tags?
	if (verbose && any(!output[[2]])) { 
        warning(paste("max iteractions exceeded for", sum(!output[[2]]), "tags", sep=" "))
    }

	output[[1]]
}
