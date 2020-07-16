maPlot <- function(x,y, logAbundance=NULL, logFC=NULL, normalize=FALSE, plot.it=TRUE, smearWidth=1, col=NULL, allCol="black", lowCol="orange", deCol="red", de.tags=NULL, smooth.scatter=FALSE, lowess=FALSE, ...)
#  Low-level function for creating an MA-plot for DGE data.
#  Created by Mark Robinson. Last modified by Davis McCarthy, 19 November 2010.
#  Edits by Gordon Smyth 20 March 2011.
{
	if( !is.null(logAbundance) && !is.null(logFC) ) {
		A <- logAbundance
		M <- logFC
		w <- v <- rep(FALSE, length(A))
		w <- A < -25+log2(1e6) # logCPM instead of logConc
		if( any(w) ) {
			shift <- max(abs(M[w])) - max(abs(M[!w]))
			A[w] <- min(A[!w]) - runif(sum(w),min=0,max=smearWidth)
			M[w] <- sign(M[w]) * (abs(M[w]) - shift)
		}
		
	} else {
		if(normalize) {
			x <- x/sum(x)
			y <- y/sum(y)
		}
		A <- (log2(x)+log2(y))/2
		M <- log2(y) - log2(x)
		w <- x==min(x) | y==min(y)
		if( any(w) ) {
			A[w] <- min(A[!w]) - runif(sum(w),min=0,max=smearWidth)
			M[w] <- log2(y[w]+min(y[!w])) - log2(x[w]+min(x[!w]))
		}
	}
	qs <- quantile(M, c(0.05,0.95))
	range <- qs[2]-qs[1]
	v <- (M < (median(M) - 5*range)) | (M > (median(M) + 5*range))
	if( any(v) ) {
		M[v] <- sign(M[v]) * (max(abs(M[!v])) + 0.5*range)
	}
	if(plot.it) {
		if( is.null(col) ) {
	  	col <- rep(allCol, length(A))
	  	if( any(w) | any(v) )
		  	col[w | v] <- lowCol
		}
		if(smooth.scatter) {
			smoothScatter(A, M, col=col, ...)
			grid()
			if( any(w) | any(v) )
				points(A[w | v], M[w | v], col=lowCol, ...)
		}
		else
			plot(A,M,col=col,...)
		points(A[de.tags],M[de.tags],col=deCol,...)
		if(lowess) {
			keep <- A > min(A[!(w |v)]) + 1
			low <- lowess(A[keep],M[keep], f=1/4)
			lines(low,col="red",lwd=4)
		}
	}
	invisible(list(A=A,M=M,w=w,v=v))
}

