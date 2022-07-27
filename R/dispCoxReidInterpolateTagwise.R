dispCoxReidInterpolateTagwise <- function(y, design, offset=NULL, dispersion, trend=TRUE, AveLogCPM=NULL, min.row.sum=5, prior.df=10, span=0.3, grid.npts=11, grid.range=c(-6,6),weights=NULL)
#	Estimate tagwise NB dispersions
#	using weighted Cox-Reid Adjusted Profile-likelihood
#	and cubic spline interpolation over a tagwise grid.
#	Yunshun Chen and Gordon Smyth
#	Created August 2010. Last modified 11 March 2013.
{
#	Check y
	y <- as.matrix(y)
	ntags <- nrow(y)
	nlibs <- ncol(y)

#	Check design
	design <- as.matrix(design)
	if(!is.fullrank(design)) stop("design matrix must be full column rank")
	ncoefs <- ncol(design)
	if(ncoefs >= nlibs) stop("no residual degrees of freedom")

#	Check offset
	lib.size <- NULL
	if(is.null(offset)) {
		lib.size <- colSums(y)
		offset <- log(lib.size)
	}
	offset <- expandAsMatrix(offset,dim(y))

#	Check AveLogCPM
	if(is.null(AveLogCPM)) AveLogCPM <- aveLogCPM(y,lib.size=lib.size)

#	Check dispersion
	ldisp <- length(dispersion)
	if(ldisp==1) {
		dispersion <- rep(dispersion,ntags)
	} else {
		if(ldisp != ntags) stop("length of dispersion doesn't match nrow(y)")
	}

#	Apply min.row.sum and use input dispersion for small count tags
	i <- (rowSums(y) >= min.row.sum)
	if(any(!i)) {
		if(any(i)) dispersion[i] <- Recall(y=y[i,,drop=FALSE],design=design,offset=offset[i,,drop=FALSE],dispersion=dispersion[i],AveLogCPM=AveLogCPM[i],grid.npts=grid.npts,min.row.sum=0,prior.df=prior.df,span=span,trend=trend,weights=weights[i,,drop=FALSE])
		return(dispersion)
	}

#	Posterior profile likelihood
	prior.n <- prior.df/(nlibs-ncoefs)
	spline.pts <- seq(from=grid.range[1],to=grid.range[2],length.out=grid.npts)
	apl <- matrix(0, nrow=ntags, ncol=grid.npts)
	for(i in 1:grid.npts){
		spline.disp <- dispersion * 2^spline.pts[i]
		apl[,i] <- adjustedProfileLik(spline.disp, y=y, design=design, offset=offset, weights=weights)
	}
	if(trend) {
		o <- order(AveLogCPM)
		oo <- order(o)
		width <- floor(span*ntags)
		apl.smooth <- movingAverageByCol(apl[o,],width=width)[oo,]
	} else {
		apl.smooth <- matrix(colMeans(apl),ntags,grid.npts,byrow=TRUE)
	}
	apl.smooth <- (apl+prior.n*apl.smooth)/(1+prior.n)

#	Tagwise maximization
	d <- maximizeInterpolant(spline.pts, apl.smooth)
	dispersion * 2^d
}

