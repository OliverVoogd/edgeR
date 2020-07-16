estimateGLMTrendedDisp <- function(y, ...) 
UseMethod("estimateGLMTrendedDisp")

estimateGLMTrendedDisp.DGEList <- function(y, design=NULL, method="auto", ...)
{
	if(is.null(y$AveLogCPM)) y$AveLogCPM <- aveLogCPM(y)

#	Check method
	method <- match.arg(method, c("auto", "bin.spline", "bin.loess", "power", "spline"))
	ntags <- nrow(y$counts)
	if(method=="auto"){
		if(ntags < 200) {
			method <- "power"
		} else {
			method <- "bin.spline"
		}
	}
	y$trend.method <- method

	d <- estimateGLMTrendedDisp(y=y$counts, design=design, offset=getOffset(y), AveLogCPM=y$AveLogCPM, method=method, weights=y$weights, ...)
	y$trended.dispersion <- d
	y
}

estimateGLMTrendedDisp.default <- function(y, design=NULL, offset=NULL, AveLogCPM=NULL, method="auto", weights=NULL, ...)
{
#	Check y
	y <- as.matrix(y)
	ntags <- nrow(y)
	if(ntags==0) return(numeric(0))
	nlibs <- ncol(y)

#	Check design
	if(is.null(design)) {
		design <- matrix(1,ncol(y),1)
		rownames(design) <- colnames(y)
		colnames(design) <- "Intercept"
	} else {
		design <- as.matrix(design)
	}
	if(ncol(design) >= ncol(y)) {
		warning("No residual df: cannot estimate dispersion")
		return(NA,ntags)
	}

#	Check offset
	if(is.null(offset)) offset <- log(colSums(y))

#	Check AveLogCPM
	if(is.null(AveLogCPM)) AveLogCPM <- aveLogCPM(y,offset=offset,weights=weights)

#	Check method	
	method <- match.arg(method,c("auto","bin.spline","bin.loess","power","spline"))
	if(method=="auto"){
		if(ntags < 200) {
			method <- "power"
		} else {
			method <- "bin.spline"
		}
	}

#	Call lower-level function
	trend <- switch(method,
		bin.spline=dispBinTrend(y, design, offset=offset, method.trend="spline", AveLogCPM=AveLogCPM, weights=weights, ...),
		bin.loess=dispBinTrend(y, design, offset=offset, method.trend="loess", AveLogCPM=AveLogCPM, weights=weights, ...),
		power=dispCoxReidPowerTrend(y, design, offset=offset, AveLogCPM=AveLogCPM, ...),
		spline=dispCoxReidSplineTrend(y, design, offset=offset, AveLogCPM=AveLogCPM, ...)
	)

	trend$dispersion
}

