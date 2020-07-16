cpmByGroup <- function(y, ...)
UseMethod("cpmByGroup")

cpmByGroup.DGEList <- function(y, group=NULL, dispersion=NULL, ...)
#	Counts per million averaged by group
#	Gordon Smyth
#	Created 10 July 2017. Last modified 4 Nov 2018.
{
	if(is.null(group)) group <- y$samples$group
	group <- as.factor(group)

	if(is.null(dispersion)) dispersion <- getDispersion(y)
	if(is.null(dispersion)) dispersion <- 0.05
	offset <- getOffset(y)

	cpmByGroup(y$counts,group=group,dispersion=dispersion,offset=offset,weights=y$weights,...)
}

cpmByGroup.SummarizedExperiment <- function(y, group=NULL, dispersion=NULL, ...)
#	Created 03 April 2020.  Last modified 03 April 2020.
{
	y <- SE2DGEList(y)
	cpmByGroup.DGEList(y, group=group, dispersion=dispersion, ...)
}

cpmByGroup.default <- function(y, group=NULL, dispersion=0.05, offset=NULL, weights=NULL, log=FALSE, prior.count=2, ...)
#	Counts per million averaged by group
#	Gordon Smyth
#	Created 10 July 2017. Last modified 4 Nov 2018.
{
	y <- as.matrix(y)

	if(is.null(group)) {
		group <- factor(rep_len(1,ncol(y)))
		levels(group) <- "AveCPM"
	}

	if(is.null(offset)) offset <- log(colSums(y))

	if(log) {
		YP <- addPriorCount(y,offset=offset,prior.count=prior.count)
		fit <- mglmOneWay(YP$y,group=group,dispersion=dispersion,offset=YP$offset,weights=weights)
		fit$coefficients / log(2) + log2(1e6)
	} else {
		fit <- mglmOneWay(y,group=group,dispersion=dispersion,offset=offset,weights=weights)
		exp(fit$coefficients) * 1e6
	}
}
