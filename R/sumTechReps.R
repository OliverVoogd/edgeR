#  sumTechReps.R

sumTechReps <- function(x,ID=colnames(x),...) UseMethod("sumTechReps")
#  24 Sept 2010

sumTechReps.default <- function(x,ID=colnames(x),...)
#	Sum over replicate columns, for matrices
#	Yifang Hu and Gordon Smyth
#	Created 14 March 2014
{
	if(is.null(x)) return(NULL)
	x <- as.matrix(x)
	if(is.null(ID)) stop("No sample IDs")
	t(rowsum(t(x),group=ID,reorder=FALSE,na.rm=FALSE))
}

sumTechReps.DGEList <- function(x,ID=colnames(x),...)
#	Sum over replicate columns, for matrices
#	Yifang Hu and Gordon Smyth
#	Created 14 March 2014. Last modified 2 Jun 2020.
{
	d <- duplicated(ID)
	if(!any(d)) return(x)
	
	x$common.dispersion <- x$trended.dispersion <- x$tagwise.dispersion <- NULL
	x$weights <- NULL
	
	y <- x[,!d]

#	Sum counts
	y$counts <- sumTechReps.default(x$counts,ID=ID,...)

#	Sum library sizes
	y$samples$lib.size <- drop(rowsum(x$samples$lib.size,group=ID,reorder=FALSE,na.rm=FALSE))

#	Average normalization factors
	y$samples$norm.factors <- drop(rowsum(x$samples$norm.factors,group=ID,reorder=FALSE,na.rm=FALSE))
	n <- rep_len(1L,nrow(x$samples))
	n <- drop(rowsum(n,group=ID,reorder=FALSE,na.rm=FALSE))
	y$samples$norm.factors <- y$samples$norm.factors/n

	rownames(y$samples) <- colnames(y$counts)
	y
}

sumTechReps.SummarizedExperiment <- function(x, ID, ...)
#	Created 03 April 2020.  Last modified 03 April 2020.
{
	x <- SE2DGEList(x)
	sumTechReps.DGEList(x, ID=ID, ...)
}
