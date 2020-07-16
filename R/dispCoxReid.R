dispCoxReid <- function(y, design=NULL, offset=NULL, weights=NULL, AveLogCPM=NULL, interval=c(0,4), tol=1e-5, min.row.sum=5, subset=10000)
#	Cox-Reid APL estimator of common dispersion
#	Gordon Smyth, Davis McCarthy
#	26 Jan 2011.  Last modified 9 Dec 2013.
{
#	Check y
	y <- as.matrix(y)

#	Check design
	if(is.null(design)) {
		design <- matrix(1,ncol(y),1)
		rownames(design) <- colnames(y)
		colnames(design) <- "Intercept"
	} else {
		design <- as.matrix(design)
	}

#	Check offset
	if(is.null(offset)) offset <- log(colSums(y))
	offset <- expandAsMatrix(offset,dim(y))
	if(min(interval)<0) stop("please give a non-negative interval for the dispersion")

#	Apply min row count
	small.row.sum <- rowSums(y)<min.row.sum
	if(any(small.row.sum)) {
		y <- y[!small.row.sum,,drop=FALSE]
		offset <- offset[!small.row.sum,,drop=FALSE]
		weights <- weights[!small.row.sum,,drop=FALSE]
		if(!is.null(AveLogCPM)) AveLogCPM <- AveLogCPM[!small.row.sum]
	}
	if(nrow(y)<1) stop("no data rows with required number of counts")

#	Subsetting
	if(!is.null(subset) && subset<=nrow(y)/2) {
		if(is.null(AveLogCPM)) AveLogCPM <- aveLogCPM(y,offset=offset,weights=weights)
		i <- systematicSubset(subset,AveLogCPM)
		y <- y[i,,drop=FALSE]
		offset <- offset[i,,drop=FALSE]
		weights <- weights[i,,drop=FALSE]
	}

#	Function for optimizing
	fun <- function(par,y,design,offset,weights) {
		sum(adjustedProfileLik(par^4,y,design,offset,weights=weights))
	}

	out <- optimize(f=fun,interval=interval^0.25,y=y,design=design,offset=offset,weights=weights,maximum=TRUE,tol=tol)
	out$maximum^4
}
