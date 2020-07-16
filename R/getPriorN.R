# ESTIMATE VALUES FOR PRIOR.N

getPriorN <- function(y, design=NULL, prior.df=20)
	## Determine the appropriate prior.n value to keep the prior degrees of freedom fixed at the given level
	## Davis McCarthy.
	## Created 29 April 2011. Last modified 29 Apr 2011.
{
	if( !is(y, "DGEList") && is.null(design) )
		stop("If y is not a DGEList object then a non-null design matrix must be provided.\n")
	nlibs <- ncol(y)
	if( is.null(design) )
		npar <- nlevels(y$samples$group)
	else
		npar <- ncol(design)
	residual.df <- nlibs - npar
	prior.n <- prior.df/residual.df
	prior.n	
}


