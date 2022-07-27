dispPearson <- function(y, design=NULL, offset=NULL, min.row.sum=5, subset=10000, AveLogCPM=NULL, tol=1e-6, trace = FALSE, initial.dispersion=0.1)
#	Pearson estimator of the common dispersion
#	using Newton iteration
#	Gordon Smyth
#	23 Aug 2012. Last modified 13 Nov 2012.
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
	if(is.null(offset)) offset <- 0
	offset <- expandAsMatrix(offset,dim(y))

#	Apply row sum filter
	small.row.sum <- which(rowSums(y)<min.row.sum)
	if(length(small.row.sum)) {
		y <- y[-small.row.sum,,drop=FALSE]
		offset <- offset[-small.row.sum,,drop=FALSE]
	}
	if(nrow(y)<1) stop("no data rows with required number of counts")

#	Apply systematic subset by AveLogCPM
	if(!is.null(subset) && subset<=nrow(y)/2) {
		if(is.null(AveLogCPM))
			AveLogCPM <- aveLogCPM(y,offset=offset)
		else {
			if(length(small.row.sum)) AveLogCPM <- AveLogCPM[-small.row.sum]
		}
		i <- systematicSubset(subset,AveLogCPM)
		y <- y[i,,drop=FALSE]
		offset <- offset[i,,drop=FALSE]
	}

#	Estimate means	using initial dispersion
	fit <- glmFit(y=y,design=design,dispersion=initial.dispersion,offset=offset,prior.count=0)
	mu <- fit$fitted.values

#	Newton iteration for dispersion, keeping means fixed
	nlibs <- ncol(y)
	df.residual <- nlibs-ncol(design)
	one <- df.residual/nlibs
	phi <- 0
	iter <- 0
	pos <- mu>0
	y <- y[pos]
	mu <- mu[pos]
	repeat {
		iter <- iter+1
		s2 <- (y-mu)^2
		Q <- mean(s2/mu/(1+phi*mu))
		dQ <- mean(s2/(1+phi*mu)^2)
		dif <- (Q-one)/dQ
		if(dif<0) break
		phi <- phi+dif
		if(trace) cat(iter,phi,Q,dQ,dif,"\n")
		if(dif < tol) break
		if(iter > 100) {
			warning("iteration limit reached")
			break
		}
	}
	phi
}
