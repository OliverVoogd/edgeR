plotMD.DGEList <- function(object, column=1, xlab="Average log CPM (this sample and others)", ylab="log-ratio (this sample vs others)", main=colnames(object)[column], status=object$genes$Status, zero.weights=FALSE, prior.count=3, ...)
#	Mean-difference plot with color coding for controls
#	Gordon Smyth
#	Created 24 June 2015. Last modified 7 Aug 2019.
{
	nlib <- ncol(object)
	if(nlib < 2L) stop("Need at least two columns")

#	Convert column to integer if not already
	j <- 1:nlib
	names(j) <- colnames(object)
	column <- j[column[1]]

	logCPM <- cpm(object, log=TRUE, prior.count=prior.count)
	AveOfOthers <- rowMeans(logCPM[,-column,drop=FALSE],na.rm=TRUE)
	Diff <- logCPM[,column]-AveOfOthers
	Mean <- (logCPM[,column]+AveOfOthers)/2

	if(!zero.weights && !is.null(object$weights)) {
		w <- as.matrix(object$weights)[,column]
		Diff[ is.na(w) | (w <= 0) ] <- NA_real_
	}

	plotWithHighlights(x=Mean,y=Diff,xlab=xlab,ylab=ylab,main=main,status=status,...)
}

plotMD.SummarizedExperiment <- function(object, column=1, xlab="Average log CPM (this sample and others)", ylab="log-ratio (this sample vs others)", zero.weights=FALSE, prior.count=3, ...)
#	Created 03 April 2020.  Last modified 03 April 2020.
{
	y <- SE2DGEList(y)
	plotMD.DGEList(y, column=column, xlab=xlab, ylab=ylab, main=colnames(y)[column], status=y$genes$Status, zero.weights=zero.weights, prior.count=prior.count, ...)
}

plotMD.DGEGLM <- function(object, column=ncol(object), coef=NULL, xlab="Average log CPM", ylab="log-fold-change", main=colnames(object)[column], status=object$genes$Status, zero.weights=FALSE, ...)
#	Mean-difference plot with color coding for controls
#	Gordon Smyth
#	Created 24 June 2015. Last modified 7 Aug 2019.
{
	if(!is.null(coef)) column <- coef
	if(is.null(object$AveLogCPM)) stop("AveLogCPM component is absent.")
	logFC <- as.matrix(object$coefficients)[,column]
	if(!zero.weights && !is.null(object$weights)) {
		w <- as.matrix(object$weights)[,column]
		logFC[ is.na(w) | (w <= 0) ] <- NA_real_
	}
	plotWithHighlights(x=object$AveLogCPM,y=logFC,xlab=xlab,ylab=ylab,main=main,status=status,...)
}

plotMD.DGELRT <- function(object, xlab="Average log CPM", ylab="log-fold-change", main=object$comparison, status=object$genes$Status, contrast=1, adjust.method="BH", p.value=0.05, ...)
#	Mean-difference plot with color coding for controls
#	Gordon Smyth
#	Created 24 June 2015. Last modified 17 June 2019.
{
	logFC <- object$table$logFC
	FTest <- is.null(logFC)
	
	if(is.null(status))
		status <- decideTestsDGE(object, adjust.method=adjust.method, p.value=p.value)

#	Multiple contrasts
	if(FTest) {
		sel <- grep("^logFC", names(object$table))[contrast]
		if(is.na(sel)) stop("Selected contrast does not exist.")
		logFC <- object$table[, sel]
		contrast.name <- gsub("logFC[.]", "", names(object$table)[sel])
		main <- paste0("Contrast ", contrast.name)
	}
	plotWithHighlights(x=object$table$logCPM,y=logFC,xlab=xlab,ylab=ylab,main=main,status=status,...)
}

plotMD.DGEExact <- function(object, xlab="Average log CPM", ylab="log-fold-change", main=NULL, status=object$genes$Status, adjust.method="BH", p.value=0.05, ...)
#	Mean-difference plot with color coding for controls
#	Gordon Smyth
#	Created 24 June 2015.  Last modified 7 Feb 2017.
{
	if(is.null(status))
		status <- decideTestsDGE(object, adjust.method=adjust.method, p.value=p.value)
	if(is.null(main)) main <- paste(object$comparison[2],"vs",object$comparison[1])
	plotWithHighlights(x=object$table$logCPM,y=object$table$logFC,xlab=xlab,ylab=ylab,main=main,status=status,...)
}
