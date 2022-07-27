dispCoxReidSplineTrend <- function(y, design, offset=NULL, df = 5, subset=10000, AveLogCPM=NULL, method.optim="Nelder-Mead", trace=0)
#	Estimate spline trend dispersion
#	Gordon Smyth, Yunshun Chen, Davis McCarthy
#	Created 16 Dec 2010.  Last modified 1 Jun 2020.
{
	y <- as.matrix(y)
	nlibs <- ncol(y)
	ntags <- nrow(y)
	lib.size <- colSums(y)
	if(is.null(offset)) offset <- 0
	offset <- expandAsMatrix(offset,dim(y))
	if(is.null(AveLogCPM)) AveLogCPM <- aveLogCPM(y,offset=offset)
	method.optim <- match.arg(method.optim, c("Nelder-Mead", "BFGS"))

	all.zero <- rowSums(y)==0
#	if(any(all.zero)) warning("Some rows of count matrix are all zero. These rows are ignored in dispersion calculation.")			   
	abundance.full <- AveLogCPM
	abundance.nonzero <- AveLogCPM[!all.zero]
	i <- systematicSubset(subset, abundance.nonzero)
	offset.nonzero <- offset[!all.zero,]
	y.nonzero <- y[!all.zero,]
	
#	Spline basis
	if(!requireNamespace("splines",quietly=TRUE)) stop("splines required but is not available")
	p1 <- (1:(df-1))/df
	knots1 <- quantile(abundance.nonzero,probs=p1)
	r <- range(abundance.nonzero)
	knots2 <- r[1]+p1*(r[2]-r[1])
	knots <- 0.3*knots1+0.7*knots2
	X <- cbind(1, splines::ns(abundance.nonzero, df, knots=knots))
	
	fun <- function(par,y,design,offset,abundance,X) {
		eta <- X %*% par
		dispersion <- exp(eta - abundance)
		tryCatch(-sum(adjustedProfileLik(as.vector(dispersion),y,design,offset)),error=function(e) 1e10)
	}

	par0 <- rep(0,df+1)
	par0[1] <- median(abundance.nonzero[i]) + log(0.1)
	out <- optim(par0,fun,y=y.nonzero[i,],design=design,offset=offset.nonzero[i,],abundance=abundance.nonzero[i],X=X[i,],control=list(trace=trace),method=method.optim)
	disp <- rep_len(NA_real_, ntags)
	disp.nonzero <-  as.vector(exp(X %*% out$par - abundance.nonzero))
	disp[all.zero] <- disp.nonzero[which.min(abundance.nonzero)]
	disp[!all.zero] <- disp.nonzero
	out$dispersion <- disp
	out$AveLogCPM <- AveLogCPM
	out
}
