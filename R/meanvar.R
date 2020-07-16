binMeanVar <- function(x, group, nbins=100L, common.dispersion=FALSE, object=NULL)
# Function to bin DGE data based on abundance and calculate the mean and pooled variance for each tag, as well as the average mean and variance for each bin. Allows us to investigate the mean-variance relationship in the data.
# Expect x to be a matrix of counts or pseudocounts---pseudocounts preferable as this adjusts for library size.
# Created by Davis McCarthy
# Created 23 Sep 2010. Last modified 7 Aug 2019.
{
	x <- as.matrix(x)
	group <- as.factor(group)
	ntags <- nrow(x)
	nlibs <- ncol(x)
	means <- rowMeans(x)

	if(nlevels(group) > 1L) {
		design <- model.matrix(~group)
		vars <- lmFit(x,design)$sigma^2
	} else {
		vars <- rowSums((x-means)^2)/(nlibs-1)
	}

	bins <- var.bins <- means.bins <- vector("list", nbins)
	o <- order(means)
	ntagsinbin <- floor(ntags / nbins)
	
	if(common.dispersion) {
		comdisp.bin <- rep_len(NA_real_, nbins)
		dispersions <- rep_len(NA_real_, nrow(x))
	} else
	dispersions <- NULL
	for(i in 1:nbins){
		if( i==nbins )
			bins[[i]] <- o[ (1 + (i-1)*ntagsinbin):ntags]
		else
			bins[[i]] <- o[ (1 + (i-1)*ntagsinbin):( i*ntagsinbin)]
		
		means.bins[[i]] <- means[bins[[i]]]
		var.bins[[i]] <- vars[bins[[i]]]
		if(common.dispersion) {
			if(!is.null(object)) {
				comdisp.bin[i] <- estimateCommonDisp(object[bins[[i]],], rowsum.filter=0)$common.dispersion
				dispersions[bins[[i]]] <- comdisp.bin[i]
			}
		}
	}
#	Take the averages of mean and variance for each bin on the square-root scale - more appropriate for the count data and reduces upward bias
	sqrt.means <- lapply(means.bins, sqrt)
	sqrt.vars <- lapply(var.bins, sqrt)
	ave.means <- sapply(sqrt.means, mean)
	ave.means <- ave.means^2
	ave.vars <- sapply(sqrt.vars, mean)
	ave.vars <- ave.vars^2
	if(common.dispersion)
		comdisp.vars <- ave.means + comdisp.bin * ave.means^2
	else {
		comdisp.vars <- NULL
		comdisp.bin <- NULL
	}
	list(avemeans=ave.means,avevars=ave.vars,bin.means=means.bins, bin.vars=var.bins, means=means, vars=vars, common.dispersion.vars=comdisp.vars, binned.common.dispersion=comdisp.bin, dispersions=dispersions, bins=bins)
}

plotMeanVar <- function(object, meanvar=NULL, show.raw.vars=FALSE, show.tagwise.vars=FALSE, show.binned.common.disp.vars=FALSE, show.ave.raw.vars=TRUE, scalar=NULL, NBline=FALSE, nbins=100, log.axes="xy", xlab=NULL, ylab=NULL, ...)
# Creates a mean-variance plot (with binned values) for a given DGEList object
# Uses the binMeanVar function
# Created by Davis McCarthy
# Created 23 Sep 2010. Last modified 7 Aug 2019.
{
	if(!is(object,"DGEList")) stop("This function requires a DGEList object")
	if(!is.null(meanvar))
		if(is.null(meanvar$means) | is.null(meanvar$vars) | is.null(meanvar$avemeans) | is.null(meanvar$avevars)) {
			message("Cannot extract all required elements of meanvar object, so will recompute it.")
			meanvar <- NULL
		}
	if(!is.null(scalar)) {
		if( length(scalar) != ncol(object$counts) )
			stop("The supplied argument scalar must have length equal to the number of columns of the count matrix in the DGEList object.")
	} else {
		offset <- getOffset(object)
		scalar <- exp(offset)/exp(mean(offset))
	}
	scalingmatrix <- expandAsMatrix(scalar, dim(object$counts))
	x <- object$counts/scalingmatrix
	if(NBline | show.tagwise.vars) {
		common.dispersion <- object$common.dispersion
		tagwise.dispersion <- object$tagwise.dispersion
		if(is.null(common.dispersion))
			stop("Could not extract common dispersion. Try running estimateCommonDisp or estimateGLMCommonDisp on the DGEList object before plotMeanVar.\n")
	}
	if(show.binned.common.disp.vars) {
		com.disp <- TRUE
		meanvar.in <- object
	} else {
		com.disp <- FALSE
		meanvar.in <- NULL
	}
	if( is.null(meanvar) ) {
		if( !is.null(object$logCPM) )
			meanvar <- binMeanVar(x, group=object$samples$group, nbins=nbins, common.dispersion=com.disp, object=meanvar.in)
		else
			meanvar <- binMeanVar(x, group=object$samples$group, nbins=nbins, common.dispersion=com.disp, object=meanvar.in)
	}
	if(show.tagwise.vars) {
		if( is.null(object$tagwise.dispersion) )
			stop("Cannot extract tagwise dispersions. Try running estimateTagwiseDisp() or estimateGLMTagwiseDisp on your object first.")
		tagvars <- meanvar$means + meanvar$means^2*tagwise.dispersion
	}
	## Averaging of means and variances in bins is now done on the square root scale to get less upward bias (done within binMeanVar) - more appropriate for count data
	avemeans <- meanvar$avemeans
	avevars <- meanvar$avevars
	
	common.dispersion.vars <- meanvar$common.dispersion.vars
	
	## Having done the necessary calculations, now do the plotting
	if(is.null(xlab))
		xlab <- "Mean gene expression level (log10 scale)"
	if(is.null(ylab))
		ylab <- "Pooled gene-level variance (log10 scale)"
	if(show.raw.vars) {
		plot(meanvar$means, meanvar$vars, log=log.axes, col="gray60", cex=0.6, xlab=xlab, ylab=ylab, plot.first=grid(), ...)
		if(show.tagwise.vars)
		  points(meanvar$means, tagvars, col="lightskyblue", cex=0.6)
		if(show.ave.raw.vars)
			points(avemeans, avevars, pch="x", col="darkred", cex=1.5)
		if(show.binned.common.disp.vars)
			points(avemeans, meanvar$common.dispersion.vars, pch="x", col="firebrick2", cex=1.5)
	}
	else {
		if(show.tagwise.vars) {
			plot(meanvar$means, tagvars, col="lightskyblue", log=log.axes, cex=0.6, xlab=xlab, ylab=ylab, plot.first=grid(), ...)
			if(show.ave.raw.vars)
				points(avemeans, avevars, pch="x", col="darkred", cex=1.5)
			if(show.binned.common.disp.vars)
				points(avemeans, meanvar$common.dispersion.vars, pch="x", col="firebrick2", cex=1.5)
		}
		else {
			if( any(!is.finite(avevars)) )
				maxy <- max(meanvar$vars)
			else
				maxy <- max(avevars)
			if(show.ave.raw.vars) {
				plot(avemeans, avevars, pch="x", col="darkred", cex=1.5, ylim=c(0.1,maxy), log=log.axes, xlab=xlab, ylab=ylab, plot.first=grid(), ...)
				if(show.binned.common.disp.vars)
					points(avemeans, meanvar$common.dispersion.vars, pch="x", col="firebrick2", cex=1.5)
			}
			else {
				plot(avemeans, meanvar$common.dispersion.vars, pch="x", col="firebrick2", cex=1.5, ylim=c(0.1,maxy), log=log.axes, xlab=xlab, ylab=ylab, plot.first=grid(), ...)
				if(show.ave.raw.vars)
					points(avemeans, avevars, pch="x", col="darkred", cex=1.5)
			}
		}
	}
	abline(0,1,lwd=2)
	if(NBline) {
		if(length(common.dispersion)==1)
			curve(x + common.dispersion*x^2, from=0.01, to=100000, col="dodgerblue3", lwd=4, add=TRUE)
		else {
			o <- order(meanvar$means)
			nb.var <- meanvar$means + (meanvar$means^2)*common.dispersion
			lines(meanvar$means[o], nb.var[o], col="dodgerblue3", lwd=4)
		}
	}
	invisible(meanvar)
}

