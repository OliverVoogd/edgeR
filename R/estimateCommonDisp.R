# Estimate common dispersion using exact conditional likelihood

estimateCommonDisp <- function(y, ...)
UseMethod("estimateCommonDisp")

estimateCommonDisp.DGEList <- function(y, tol=1e-06, rowsum.filter=5, verbose=FALSE, ...)
# Yunshun Chen. Created 7 Aug 2019.
{
	y <- validDGEList(y)
	group <- y$samples$group
	lib.size <- y$samples$lib.size * y$samples$norm.factors
	
	if( all(tabulate(group)<=1) ) {
		warning("There is no replication, setting dispersion to NA.")
		y$common.dispersion <- NA_real_
		return(y)
	}

	out <- estimateCommonDisp(y$counts, group=group, lib.size=lib.size, tol=tol, rowsum.filter=rowsum.filter, verbose=verbose, ...)	
	y$common.dispersion <- out
	y <- equalizeLibSizes(y, dispersion=out)
	y$AveLogCPM <- aveLogCPM(y, dispersion=out)
	y
}


estimateCommonDisp.default <- function(y, group=NULL, lib.size=NULL, tol=1e-06, rowsum.filter=5, verbose=FALSE, ...)
#	Davis McCarthy, Mark Robinson, Gordon Smyth.
#	Created 2009. Last modified 18 March 2016.
{
#	Check y
	y <- as.matrix(y)
	ntags <- nrow(y)
	if(ntags==0L) stop("No data rows")
	nlibs <- ncol(y)
	if(nlibs < 2L) stop("Need at least two libraries")

#	Check group
	if(is.null(group)) group <- rep(1, nlibs)
	if(length(group)!=nlibs) stop("Incorrect length of group.")
	group <- dropEmptyLevels(group)

#	Check lib.size
	if(is.null(lib.size)) {
		lib.size <- colSums(y)
	} else {
		if(length(lib.size)!=nlibs) stop("Incorrect length of lib.size.")
	}

#	Filter low count genes
	sel <- rowSums(y) > rowsum.filter
	if(!sum(sel)) stop("No genes satisfy rowsum filter")
	
#	Start from small dispersion
	disp <- 0.01
	for(i in 1:2) {
		out <- equalizeLibSizes(y, group=group, dispersion=disp, lib.size=lib.size)
		y.pseudo <- out$pseudo.counts[sel, , drop=FALSE]
		y.split <- splitIntoGroups(y.pseudo, group=group)
		delta <- optimize(commonCondLogLikDerDelta, interval=c(1e-4,100/(100+1)), tol=tol, maximum=TRUE, y=y.split, der=0)
		delta <- delta$maximum
		disp <- delta/(1-delta)
	}
	if(verbose) cat("Disp =",round(disp,5),", BCV =",round(sqrt(disp),4),"\n")
	
	common.dispersion <- disp
	common.dispersion
}

