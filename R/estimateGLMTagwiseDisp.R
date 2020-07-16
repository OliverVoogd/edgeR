estimateGLMTagwiseDisp <- function(y, ...) 
UseMethod("estimateGLMTagwiseDisp")

estimateGLMTagwiseDisp.DGEList <- function(y, design=NULL, prior.df=10, trend=!is.null(y$trended.dispersion), span=NULL, ...)
{
#	Find appropriate dispersion
	if(trend) {
		dispersion <- y$trended.dispersion
		if(is.null(dispersion)) stop("No trended.dispersion found in data object. Run estimateGLMTrendedDisp first.")
	} else {
		dispersion <- y$common.dispersion
		if(is.null(dispersion)) stop("No common.dispersion found in data object. Run estimateGLMCommonDisp first.")
	}

	if(is.null(y$AveLogCPM)) y$AveLogCPM <- aveLogCPM(y)
	ntags <- nrow(y$counts)
	if(is.null(span)) if(ntags>10) span <- (10/ntags)^0.23 else span <- 1
	y$span <- span
	
	d <- estimateGLMTagwiseDisp(y=y$counts, design=design, offset=getOffset(y), dispersion=dispersion, trend=trend, span=span, prior.df=prior.df, AveLogCPM=y$AveLogCPM, weights=y$weights, ...)
	y$prior.df <- prior.df
	y$tagwise.dispersion <- d
	y
}

estimateGLMTagwiseDisp.default <- function(y, design=NULL, offset=NULL, dispersion, prior.df=10, trend=TRUE, span=NULL, AveLogCPM=NULL, weights=NULL, ...)
#	Created 24 March 2011. Last modified 7 Aug 2019.
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
		warning("No residual df: setting dispersion to NA")
		return(rep_len(NA_real_,ntags))
	}

#	Check offset
	if(is.null(offset)) offset <- log(colSums(y))

#	Check span
#	span can be chosen smaller when ntags is large
	if(is.null(span)) if(ntags>10) span <- (10/ntags)^0.23 else span <- 1

#	Check AveLogCPM
	if(is.null(AveLogCPM)) AveLogCPM <- aveLogCPM(y, offset=offset, weights=weights)

#	Call Cox-Reid grid method
	tagwise.dispersion <- dispCoxReidInterpolateTagwise(y, design, offset=offset, dispersion, trend=trend, prior.df=prior.df, span=span, AveLogCPM=AveLogCPM, weights=weights,...)
	tagwise.dispersion
}
