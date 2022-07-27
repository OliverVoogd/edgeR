rowsum.DGEList <- function (x, group, reorder=FALSE, na.rm=FALSE, ...)
#	Sum counts by groups of rows/genes and
#	return a DGEList with a row for each level of 'group'.
#	Gordon Smyth
#	Created 22 Feb 2018. Last modified 23 May 2018.
{
	isdupgrp <- duplicated(group)
	x2 <- x[!isdupgrp,]
	x2$counts <- rowsum(x$counts,group=group,reorder=FALSE,na.rm=na.rm,...)
	if(!is.null(x$genes)) {
#		Keep those columns of x$genes that contain group-level annotation
		no <- logical(nrow(x))
		isdupall <- vapply(x$genes,duplicated,no)[isdupgrp,,drop=FALSE]
		isgenelevel <- (colSums(isdupall) == nrow(isdupall))
		x2$genes <- x2$genes[,isgenelevel,drop=FALSE]
		row.names(x2$genes) <- row.names(x2$counts)
	}
	if(reorder) {
		o <- order(row.names(x2))
		x2 <- x2[o,]
	}
	x2
}

rowsum.SummarizedExperiment <- function(x, group, reorder=FALSE, na.rm=FALSE, ...)
#	Created 03 April 2020.  Last modified 03 April 2020.
{
	x <- SE2DGEList(x)
	rowsum.DGEList(x, group=group, reorder=reorder, na.rm=na.rm, ...)
}
