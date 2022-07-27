predFC <- function(y,design,prior.count=0.125,offset=NULL,dispersion=NULL,weights=NULL,...) 
UseMethod("predFC")

predFC.DGEList <- function(y,design,prior.count=0.125,offset=NULL,dispersion=NULL,weights=NULL,...)
{
	if(is.null(offset)) offset <- getOffset(y)
	if(is.null(dispersion)) dispersion <- getDispersion(y)
	if(is.null(dispersion)) {
		dispersion <- 0
		message("dispersion set to zero")
	}
	predFC.default(y=y$counts,design=design,prior.count=prior.count,offset=offset,dispersion=dispersion,weights=weights,...)
}

predFC.SummarizedExperiment <- function(y,design,prior.count=0.125,offset=NULL,dispersion=NULL,weights=NULL,...)
#	Created 03 April 2020.  Last modified 03 April 2020.
{
	y <- SE2DGEList(y)
	predFC.DGEList(y, design=design, prior.count=prior.count, offset=offset, dispersion=dispersion, weights=weights, ...)
}

predFC.default <- function(y,design,prior.count=0.125,offset=NULL,dispersion=0,weights=NULL,...)
#	Shrink log-fold-changes towards zero by augmenting data counts
#	Gordon Smyth and Belinda Phipson
#	17 Aug 2011.  Last modified 9 July 2017.
{
#	Add prior counts in proportion to library sizes
	out <- addPriorCount(y, offset=offset, prior.count=prior.count)

#	Check design
	design <- as.matrix(design)

#	Return matrix of coefficients on log2 scale
	g <- glmFit(out$y,design,offset=out$offset,dispersion=dispersion,prior.count=0,weights=weights,...)
	g$coefficients/log(2)
}

