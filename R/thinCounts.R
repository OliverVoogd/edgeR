thinCounts <- function(x,prob=NULL,target.size=min(colSums(x)))
#	Binomial or multinomial thinning of counts
#	Gordon Smyth
#	23 March 2011.  Last revised 23 Nov 2011.
{
	if(!is.null(prob)) {
		x[] <- rbinom(length(x),size=x,prob=prob)
	} else {
		x <- as.matrix(x)
		target.size <- rep.int(target.size,ncol(x))
		actual.size <- colSums(x)
		if(any(target.size>actual.size)) stop("target.size bigger than actual size")
		for (j in 1:ncol(x)) {
			diff.size <- actual.size[j]-target.size[j]
			if(diff.size>0) x[,j] <- x[,j]-rmultinom(1,size=diff.size,prob=x[,j])
		}
		if(any(x<0)) {
			x <- pmax(x,0)
			x <- Recall(x,target.size=target.size)
		}
	}
	x
}

