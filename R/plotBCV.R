plotBCV <- function(y, xlab="Average log CPM", ylab="Biological coefficient of variation", pch=16, cex=0.2, col.common="red", col.trend="blue", col.tagwise="black", ...)
#	Plot biological coefficient of variation against average log CPM
#	Davis McCarthy, Yunshun Chen, Gordon Smyth.
#	Created 18 January 2012.  Last modified 11 March 2013.
{
#	Check y
	if(!is(y,"DGEList")) stop("y must be a DGEList.")

#	Compute AveLogCPM if not found in y
	A <- y$AveLogCPM
	if(is.null(A)) A <- aveLogCPM(y$counts, offset=getOffset(y))

#	Points to determine y axis limits
	disp <- getDispersion(y)
	if(is.null(disp)) stop("No dispersions to plot")
	if(attr(disp,"type")=="common") disp <- rep_len(disp, length(A))

#	Make plot
	plot(A, sqrt(disp), xlab=xlab, ylab=ylab, type="n", ...)
	labels <- cols <- lty <- pt <- NULL
	if(!is.null(y$tagwise.dispersion)) {
		points(A, sqrt(y$tagwise.dispersion), pch=pch, cex=cex, col=col.tagwise)
		labels <- c(labels, "Tagwise")
		cols <- c(cols, col.tagwise)
		lty <- c(lty, -1)
		pt <- c(pt, pch)
	}
	if(!is.null(y$common.dispersion)) {
		abline(h=sqrt(y$common.dispersion), col=col.common, lwd=2)
		labels <- c(labels, "Common")
		cols <- c(cols, col.common)
		lty <- c(lty, 1)
		pt <- c(pt, -1)		
	}
	if(!is.null(y$trended.dispersion)) {
		o <- order(A)
		lines(A[o], sqrt(y$trended.dispersion)[o], col=col.trend, lwd=2)
		labels <- c(labels, "Trend")
		cols <- c(cols, col.trend)
		lty <- c(lty, 1)
		pt <- c(pt, -1)		
	}
	legend("topright", legend=labels, lty=lty, pch=pt, pt.cex=cex, lwd=2, col=cols)

	invisible()
}
