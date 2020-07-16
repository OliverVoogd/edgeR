plotMDS.DGEList <- function (x,top=500,labels=NULL,pch=NULL,cex=1,dim.plot=c(1,2),ndim=max(dim.plot),gene.selection="pairwise",xlab=NULL,ylab=NULL,method="logFC",prior.count=2,plot=TRUE,...)
#	Multidimensional scaling plot of digital gene expression profiles
#	Yunshun Chen, Mark Robinson and Gordon Smyth
#	23 May 2011.  Last modified 26 Nov 2016.
{
	method <- match.arg(method, c("logfc","logFC","bcv","BCV"))
	if(method=="logfc") method <- "logFC"
	if(method=="BCV") method <- "bcv"

#	Default method is to convert to moderated logCPM and call limma plotMDS
	if(method=="logFC") {
		y <- cpm(x,log=TRUE,prior.count=prior.count)
		return(plotMDS(y,top=top,labels=labels,pch=pch,cex=cex,dim.plot=dim.plot,ndim=ndim,gene.selection=gene.selection,xlab=xlab,ylab=ylab,plot=plot,...))
	}

#	From here method="bcv"

#	Check x
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
	if(ndim < 2) stop("Need at least two dim.plot")
	if(nsamples < ndim) stop("Too few samples")
	if(nprobes < ndim) stop("Too few rows")

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

#	Multidim scaling
	a1 <- cmdscale(as.dist(dd), k = ndim)

#	Check whether dimensions have been removed (because of negative eigenvalues)
#	Add random variate if necessary
	ndiff <- ndim-ncol(a1)
	if(ndiff > 0) a1 <- cbind(a1, matrix(runif(ndiff*nsamples, -1e-6, 1e-6), ncol=ndiff, nrow=nsamples))

	mds <- new("MDS",list(dim.plot=dim.plot,distance.matrix=dd,cmdscale.out=a1,top=top))
	mds$x <- a1[,dim.plot[1]]
	mds$y <- a1[,dim.plot[2]]
	mds$axislabel <- "BCV distance"
	if(plot)
		plotMDS(mds,labels=labels,pch=pch,cex=cex,xlab=xlab,ylab=ylab,...)
	else
		mds
}

plotMDS.SummarizedExperiment <- function(x, top=500, labels=NULL, pch=NULL, cex=1, dim.plot=c(1,2), ndim=max(dim.plot), gene.selection="pairwise", xlab=NULL, ylab=NULL, method="logFC", prior.count=2, plot=TRUE, ...)
#	Created 03 April 2020.  Last modified 03 April 2020.
{
	x <- SE2DGEList(x)
	plotMDS.DGEList(x, top=top, labels=labels, pch=pch, cex=cex, dim.plot=dim.plot, ndim=ndim, gene.selection=gene.selection, xlab=xlab, ylab=ylab, method=method, prior.count=prior.count, plot=plot, ...)
}
