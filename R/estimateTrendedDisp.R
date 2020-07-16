# Estimate trended dispersions using exact conditional likelihood

estimateTrendedDisp <- function(y, ...)
UseMethod("estimateTrendedDisp")

estimateTrendedDisp.DGEList <- function(y, method="bin.spline", df=5, span=2/3, ...)
# Yunshun Chen. Created 18 March 2016.
{
	y <- validDGEList(y)
	group <- y$samples$group
	lib.size <- y$samples$lib.size * y$samples$norm.factors
	if(is.null(y$AveLogCPM)) y$AveLogCPM <- aveLogCPM(y)
	
	out <- estimateTrendedDisp(y$counts, group=group, lib.size=lib.size, AveLogCPM=y$AveLogCPM, method=method, df=df, span=span)
	y$trended.method <- method
	y$trended.dispersion <- out
	y
}


estimateTrendedDisp.default <- function(y, group=NULL, lib.size=NULL, AveLogCPM=NULL, method="bin.spline", df=5, span=2/3, ...)
# Yunshun Chen, Gordon Smyth.
# Created 2 Feb 2012, last modified on 16 May 2019.
{
#	Check y
	y <- as.matrix(y)
	ntags <- nrow(y)
	nlibs <- ncol(y)

#	Check group
	if(is.null(group)) group <- rep(1, nlibs)
	if(length(group)!=nlibs) stop("Incorrect length of group.")
	group <- dropEmptyLevels(group)

#	Check lib.size
	if(is.null(lib.size)) lib.size <- colSums(y)
	if(length(lib.size)!=nlibs) stop("Incorrect length of lib.size.")

#	Check AveLogCPM
	if(is.null(AveLogCPM)) AveLogCPM <- aveLogCPM(y, lib.size=lib.size)

#	Check method
	method <- match.arg(method, c("bin.spline", "bin.loess"))
	
	nbins <- 50
	if(nbins > ntags) stop("nbins greater than number of rows of data")
	bins <- cutWithMinN(AveLogCPM, intervals=nbins, min.n=floor(ntags/nbins))
	disp.bins <- AveLogCPM.bins <- rep_len(NA_real_, nbins)
	
	for(i in 1:nbins) {
		tagsinbin <- bins$group==i
		disp.bins[i] <- estimateCommonDisp(y[tagsinbin,], group=group, lib.size=lib.size, rowsum.filter=0)
		AveLogCPM.bins[i] <- mean(AveLogCPM[tagsinbin])
	}

	if( method=="bin.spline" ) {
		if(!requireNamespace("splines",quietly=TRUE)) stop("splines required but is not available")
		p1 <- (1:(df-1))/df
		knots1 <- quantile(AveLogCPM.bins, probs=p1)
		r <- range(AveLogCPM.bins)
		knots2 <- r[1] + p1*(r[2]-r[1])
		knots <- 0.3*knots1 + 0.7*knots2
		basisbins <- splines::ns(AveLogCPM.bins, df=df, knots=knots, intercept=TRUE)
		beta <- coefficients(lm.fit(basisbins, sqrt(disp.bins)))
		basisall <- predict(basisbins, newx=AveLogCPM)
		trended.dispersion <- drop(basisall %*% beta)^2
	}
	
	if( method=="bin.loess" ) {
		fit <- loessFit(sqrt(disp.bins), AveLogCPM.bins, span=span, iterations=1)
		f <- approxfun(AveLogCPM.bins, fit$fitted, rule=2, ties=list("ordered",mean))
		trended.dispersion <- f(AveLogCPM)^2
	}
	
	trended.dispersion
}

