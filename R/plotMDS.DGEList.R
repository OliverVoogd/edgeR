plotMDS.DGEList <- function (x,top=500,labels=NULL,pch=NULL,cex=1,dim.plot=c(1,2),gene.selection="pairwise",xlab=NULL,ylab=NULL,method="logFC",prior.count=2,plot=TRUE,var.explained=TRUE,...)
#	Multidimensional scaling plot of digital gene expression profiles
#	Yunshun Chen, Mark Robinson and Gordon Smyth
#	23 May 2011.  Last modified 13 May 2021.
{
	method <- match.arg(method, c("logfc","logFC","bcv","BCV"))
	if(method=="logfc") method <- "logFC"
	if(method=="BCV") method <- "bcv"

#	Default method is to convert to moderated logCPM and call limma plotMDS
	if(method=="logFC") {
		y <- cpm(x,log=TRUE,prior.count=prior.count)
		return(plotMDS(y,top=top,labels=labels,pch=pch,cex=cex,dim.plot=dim.plot,gene.selection=gene.selection,xlab=xlab,ylab=ylab,plot=plot,var.explained=var.explained,...))
	}

#	From here method="bcv"

#	Check x
	message("Note: the bcv method is now scheduled to be removed in a future release of edgeR.")
	x$counts <- as.matrix(x$counts)
	if(!all(is.finite(x$counts))) stop("Missing or infinite counts not allowed")

	nprobes <- nrow(x)
	nsamples <- ncol(x)
	if(nsamples < 3) stop("Need at least 3 columns of data")

#	Check labels and pch
	if(is.null(pch) & is.null(labels)) {
		labels <- colnames(x)
		if(is.null(labels)) labels <- 1:nsamples
	}
	if(!is.null(labels)) labels <- as.character(labels)

#	Check dim
	ndim <- max(dim.plot)
	if(ndim < 2) stop("Need at least two dim.plot")
	if(nsamples < ndim) stop("ndim is greater than number of samples")
	if(nprobes < ndim) stop("ndim is greater than number of rows of data")

	x$samples$group <- factor(rep.int(1,nsamples))

	cn <- colnames(x)
	dd <- matrix(0,nrow=nsamples,ncol=nsamples,dimnames=list(cn,cn))	

#	Check value for top
	if(top < nprobes) { 
		twd <- estimateTagwiseDisp(estimateCommonDisp(x), grid.length = 100) 
		o <- order(twd$tagwise.dispersion, decreasing = TRUE)[1:top]
		subdata <- x$counts[o,,drop=FALSE]
	} else {
		subdata <- x$counts
	}

#	Compute pairwise BCV values
	lib.size <- x$samples$lib.size * x$samples$norm.factors
	myFun <- function(delta, y, ...) sum(condLogLikDerDelta(y, delta, ...))
	for (i in 2:(nsamples)) {
		for (j in 1:(i - 1))  {
			mm <- subdata[,c(i,j)]
			rs5 <- rowSums(mm) > 5
			lib <- lib.size[c(i, j)]
			norm <- t(t(mm)/lib) * exp(mean(log(lib)))
			delta <- optimize(myFun, interval = c(0.0001,.99), tol = 0.000001, maximum = TRUE, y = norm[rs5,], der = 0)
			dd[i, j] = sqrt( delta$maximum / (1-delta$maximum) )
		}
	}

#	Multi-dimensional scaling
	dd <- dd + t(dd)
	rm <- rowMeans(dd)
	dd <- dd - rm
	dd <- t(dd) - (rm - mean(rm))
	mds <- eigen(-dd/2, symmetric=TRUE)
	names(mds) <- c("eigen.values","eigen.vectors")

#	Make MDS object
	lambda <- pmax(mds$eigen.values,0)
	mds$var.explained <- lambda / sum(lambda)
	mds$dim.plot=dim.plot
	mds$distance.matrix.squared=dd
	mds$top=top
	mds$gene.selection=gene.selection
	mds$axislabel <- "BCV distance"
	mds <- new("MDS",unclass(mds))

#	Add coordinates for plot
	i <- dim.plot[1]
	mds$x <- mds$eigen.vectors[,i] * sqrt(lambda[i])
	if(lambda[i] < 1e-13) warning("dimension ", i, " is degenerate or all zero")
	i <- dim.plot[2]
	mds$y <- mds$eigen.vectors[,i] * sqrt(lambda[i])
	if(lambda[i] < 1e-13) warning("dimension ", i, " is degenerate or all zero")

	if(plot)
		plotMDS(mds,labels=labels,pch=pch,cex=cex,xlab=xlab,ylab=ylab,var.explained=var.explained,...)
	else
		mds
}

plotMDS.SummarizedExperiment <- function(x, top=500, labels=NULL, pch=NULL, cex=1, dim.plot=c(1,2), gene.selection="pairwise", xlab=NULL, ylab=NULL, method="logFC", prior.count=2, plot=TRUE, var.explained=TRUE, ...)
#	Created 03 April 2020.  Last modified 13 May 2021.
{
	x <- SE2DGEList(x)
	plotMDS.DGEList(x, top=top, labels=labels, pch=pch, cex=cex, dim.plot=dim.plot, gene.selection=gene.selection, xlab=xlab, ylab=ylab, method=method, prior.count=prior.count, plot=plot, var.explained=var.explained, ...)
}
