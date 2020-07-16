rpkmByGroup <- function(y, ...)
UseMethod("rpkmByGroup")

rpkmByGroup.DGEList <- function(y, group=NULL, gene.length=NULL, dispersion=NULL, ...)
#	RPKM or FPKM averaged by group
#	Gordon Smyth
#	Created 10 July 2017. Last modified 4 Nov 2018.
{
	if(is.null(group)) group <- y$samples$group
	group <- as.factor(group)

	if(is.null(gene.length)) {
#		Look for a column name containing "length"
		gene.length <- y$genes$Length
		if(is.null(gene.length)) gene.length <- y$genes$length
		if(is.null(gene.length)) {
			j <- grep("length",tolower(names(y$genes)))
			if(length(j)==1L)
				gene.length <- y$genes[,j]
			else
				stop("Gene lengths not found")
		}
	} else {
		if(is.character(gene.length)) {
#			A column name of y$genes has been specified
			gene.length <- y$genes[[gene.length[1]]]
			if(is.null(gene.length)) stop("gene.length column not found")
		}
		gene.length <- as.numeric(gene.length)
		if(length(gene.length) != nrow(y)) stop("length of gene.length doesn't match nrows of y")
	}

	if(is.null(dispersion)) dispersion <- getDispersion(y)
	if(is.null(dispersion)) dispersion <- 0.05
	offset <- getOffset(y)

	rpkmByGroup(y$counts,group=group,gene.length=gene.length,dispersion=dispersion,offset=offset,weights=y$weights,...)
}

rpkmByGroup.SummarizedExperiment <- function(y, group=NULL, gene.length=NULL, dispersion=NULL, ...)
#	Created 03 April 2020.  Last modified 03 April 2020.
{
	y <- SE2DGEList(y)
	rpkmByGroup.DGEList(y, group=group, gene.length=gene.length, dispersion=dispersion, ...)
}

rpkmByGroup.default <- function(y, group=NULL, gene.length, dispersion=0.05, offset=NULL, weights=NULL, log=FALSE, prior.count=2, ...)
#	RPKM or FPKM averaged by group
#	Gordon Smyth
#	Created 10 July 2017. Last modified 21 June 2019.
{
	if(is.null(group)) {
		group <- factor(rep_len(1,ncol(y)))
		levels(group) <- "AveRPKM"
	}

	z <- cpmByGroup(y=y,group=group,dispersion=dispersion,offset=offset,weights=weights,log=log,prior.count=prior.count, ...)
	if(log) {
		z - log2(gene.length / 1e3)
	} else {
		z / (gene.length / 1e3)
	}

}
