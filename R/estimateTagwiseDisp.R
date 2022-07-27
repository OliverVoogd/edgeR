#  Tagwise dispersion using weighted conditional likelihood empirical Bayes.

estimateTagwiseDisp <- function(y, ...) 
UseMethod("estimateTagwiseDisp")

estimateTagwiseDisp.DGEList <- function(y, prior.df=10, trend="movingave", span=NULL, method="grid", grid.length=11, grid.range=c(-6,6), tol=1e-06, verbose=FALSE, ...)
# Yunshun Chen. Created 18 March 2016.
{
	y <- validDGEList(y)
	group <- y$samples$group
	lib.size <- y$samples$lib.size * y$samples$norm.factors
	if(is.null(y$AveLogCPM)) y$AveLogCPM <- aveLogCPM(y)
	dispersion <- y$common.dispersion
	if(is.null(dispersion)) stop("No common.dispersion found in the DGEList object. Run estimateCommonDisp first.")

	y$prior.df <- prior.df
	nlibs <- ncol(y$counts)
	ngroups <- length(unique(group))
	y$prior.n <- prior.df/(nlibs - ngroups)
	if(is.null(span)) if(trend=="movingave") span <- 0.3 else span <- 0.5
	
	out <- estimateTagwiseDisp(y$counts, group=group, lib.size=lib.size, dispersion=dispersion, AveLogCPM=y$AveLogCPM, prior.df=prior.df, trend=trend, span=span, method=method, grid.length=grid.length, grid.range=grid.range, tol=tol, verbose=verbose)
	y$tagwise.dispersion <- out
	y$span <- span
	y
}

estimateTagwiseDisp.default <- function(y, group=NULL, lib.size=NULL, dispersion, AveLogCPM=NULL, prior.df=10, trend="movingave", span=NULL, method="grid", grid.length=11, grid.range=c(-6,6), tol=1e-06, verbose=FALSE, ...)
#  Davis McCarthy, Mark Robinson, Yunshun Chen, Gordon Smyth.
#  Created 2009. Last modified 1 June 2016.

#  Notes 3 July 2012:
#  - interpolating derivatives would be better than interpolating loglik values.
#  - share code with estimateGLMTagwiseDisp?
{
#	Check y
	y <- as.matrix(y)
	ntags <- nrow(y)
	nlibs <- ncol(y)

#	Check group
	if(is.null(group)) group <- rep(1, nlibs)
	if(length(group)!=nlibs) stop("Incorrect length of group.")
	group <- dropEmptyLevels(group)

#	Check lib.size
	if(is.null(lib.size)) lib.size <- colSums(y)
	if(length(lib.size)!=nlibs) stop("Incorrect length of lib.size.")

#	Check AveLogCPM
	if(is.null(AveLogCPM)) AveLogCPM <- aveLogCPM(y, lib.size=lib.size, dispersion=dispersion)

	trend <- match.arg(trend, c("none", "loess", "movingave", "tricube"))
	if(trend=="tricube") trend <- "loess"
	method <- match.arg(method, c("grid", "optimize"))

	eq <- equalizeLibSizes(y, group=group, dispersion=dispersion, lib.size=lib.size)	
	u <- splitIntoGroups(eq$pseudo.counts, group=group)
	delta <- rep(0, ntags)
	ngroups <- length(unique(group))
	prior.n <- prior.df/(nlibs - ngroups)

	if(method=="grid"){  # do spline interpolation
		if(verbose) message("Using interpolation to estimate tagwise dispersion. ")
		spline.pts <- seq(from=grid.range[1], to=grid.range[2], length.out=grid.length)
		spline.disp <- dispersion * 2^spline.pts
		grid.vals <- spline.disp/(1+spline.disp)
	
		l0 <- matrix(0, ntags, grid.length)
		for(j in 1:grid.length) for(i in 1:length(u)) l0[,j] <- condLogLikDerDelta(u[[i]], grid.vals[j], der=0L) + l0[,j]

		if(is.null(span)) if(trend=="movingave") span <- 0.3 else span <- 0.5
		m0 <- switch(trend,
 			"movingave" = {
 				o <- order(AveLogCPM)
 				oo <- order(o)
 				movingAverageByCol(l0[o,], width=floor(span*ntags))[oo,]
 			},
			"loess" = loessByCol(l0, AveLogCPM, span=span)$fitted.values,
			"none" = matrix(colMeans(l0), ntags, grid.length, byrow=TRUE)
		)
		l0a <- l0 + prior.n*m0
		d <- maximizeInterpolant(spline.pts, l0a)
		tagwise.dispersion <- dispersion * 2^d
	} else {	
		if(trend != "none") stop("optimize method doesn't allow for abundance-dispersion trend")
		if(verbose) message("Tagwise dispersion optimization begun, may be slow, progress reported every 100 tags")
		for(tag in seq_len(ntags)) {
			delta.this <- optimize(weightedCondLogLikDerDelta, interval=c(1e-4,100/(100+1)), tol=tol, maximum=TRUE, y=u, tag=tag, ntags=ntags, prior.n=prior.n, der=0L)
			delta[tag] <- delta.this$maximum
			if(verbose) if(tag%%100==0) message("tag ",tag)
		}
		tagwise.dispersion <- delta/(1-delta)
	}
	if(verbose) cat("\n")

	tagwise.dispersion
}

