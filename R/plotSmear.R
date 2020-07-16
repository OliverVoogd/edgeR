plotSmear <- function (object, pair=NULL, de.tags=NULL, xlab="Average logCPM", ylab="logFC", pch=19, cex=.2, smearWidth=.5, panel.first=grid(), smooth.scatter=FALSE, lowess=FALSE, ...)
# User-level function for creating an MA-plot for DGE data.
# Created 3 Oct 2009. Last modified 1 Jun 2020.
 
{
	if ( !inherits(object,c("DGEList", "DGELRT", "DGEExact")) ) 
		stop("Currently only supports DGEList/DGELRT/DGEExact objects as the object argument.")
	if( is(object, "DGEList") && is.null(object$samples$group) )
		stop("Cannot produce a smear plot if no experimental groups are defined. Here, object$samples$group is NULL.\n")
	if( is(object, "DGEList") ) {
		levs.group <- levels(object$samples$group)
		if(length(levs.group)==1)
			stop("Cannot produce an MA-plot with only one group. The one group defined is: ",levs.group)
		if (is.null(pair))
			pair <- levs.group[1:2]
		if( !all(pair %in% levs.group) )
			stop("At least one element of given pair is not a group.\n Groups are: ", paste(levs.group, collapse=" "), "\n")
		stopifnot(length(pair)==2)
	   
		cols1 <- pair[1]==object$samples$group
		cols2 <- pair[2]==object$samples$group
		
		lib.size <- object$samples$lib.size*object$samples$norm.factors
		x <- 1e6 * rowMeans( object$counts[,cols1,drop=FALSE] / expandAsMatrix( lib.size[cols1], dim(object$counts[,cols1,drop=FALSE])) )
		y <- 1e6 * rowMeans( object$counts[,cols2,drop=FALSE] / expandAsMatrix( lib.size[cols2], dim(object$counts[,cols2,drop=FALSE])) )
		
		ylab <- paste(ylab, ":", paste(pair[c(2,1)], collapse="-"), paste="")
		i <- match(de.tags,rownames(object$counts))
		i <- i[!is.na(i)]
		
		maPlot( x, y, xlab=xlab, ylab=ylab, pch=pch, cex=cex, smearWidth=smearWidth, de.tags=i, panel.first=panel.first, smooth.scatter=smooth.scatter, lowess=lowess, ...)
	} else {
		if(is.null(object$table$logFC))
			stop("table$logFC slot in DGELRT object is NULL. We cannot produce an MA (smear) plot if more than one coefficient from the GLM is being tested in the likelihood ratio test as this results in more one logFC value per gene---one for each coefficient.\n")
		i <- match(de.tags,rownames(object$table))
		i <- i[!is.na(i)]
		maPlot( x=NULL, y=NULL, logAbundance=object$table$logCPM, logFC=object$table$logFC, xlab=xlab, ylab=ylab, pch=pch, cex=cex, smearWidth=smearWidth, de.tags=i, panel.first=panel.first, smooth.scatter=smooth.scatter, lowess=lowess, ...)
	}
}

