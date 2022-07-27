validDGEList <- function(y)
#	Check for standard components of DGEList object
#	Gordon Smyth
#	20 Nov 2013
{
	if(is.null(y$counts)) stop("No count matrix")
	y$counts <- as.matrix(y$counts)
	nlib <- ncol(y$counts)
	if(is.null(y$samples$group)) y$samples$group <- gl(1,nlib)
	if(is.null(y$samples$lib.size)) y$samples$lib.size <- colSums(y$counts)
	if(is.null(y$samples$norm.factors)) y$samples$norm.factors <- rep.int(1,nlib)
	y
}
