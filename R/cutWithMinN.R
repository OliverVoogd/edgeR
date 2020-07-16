cutWithMinN <- function(x, intervals=2, min.n=1)
#	Cut numeric x into intervals, as equally spaced as possible subject
#	to including a minimum number of values in each interval
#	Gordon Smyth
#	7 May 2011.  Last modified 9 Jun 2020.
{
#	Check input
	x <- as.numeric(x)
	isna <- is.na(x)
	if(any(isna)) {
		group <- rep_len(NA_real_,length(x))
		out <- Recall(x=x[!isna],intervals=intervals,min.n=min.n)
		group[!isna] <- out$group
		out$group <- group
		return(out)
	}
	intervals <- as.integer(intervals)
	min.n <- as.integer(min.n)
	nx <- length(x)
	if(nx < intervals*min.n) stop("too few observations: length(x) < intervals*min.n")

	if(intervals==1) return(list(group=rep(1,nx),breaks=NA))

#	Add jittering to ensure all x are unique
#	x <- x+(1e-10)*(1:nx)/nx
	x <- x+(1e-10)*(runif(nx)-0.5)

#	Breaks equally spaced by x
	breaks.eqx <- seq(from=min(x),to=max(x),length.out=intervals+1L)
	breaks.eqx[1] <- breaks.eqx[1]-1
	breaks.eqx[intervals+1L] <- breaks.eqx[intervals+1L]+1

#	Breaks equally spaced by quantiles
	breaks.eqn <- quantile(x,probs=seq(from=0,to=1,length.out=intervals+1L))
	breaks.eqn[1] <- breaks.eqn[1]-1
	breaks.eqn[intervals+1L] <- breaks.eqn[intervals+1L]+1

#	First try equally spaced by x
	z <- cut(x,breaks=intervals,labels=FALSE)
	n <- tabulate(z)
	if(all(n>=min.n)) return(list(group=z,breaks=breaks.eqx))

#	Step down gradually
	for (i in 1:9) {
		breaks <- (i*breaks.eqn+(10-i)*breaks.eqx)/10
		z <- cut(x,breaks=breaks,labels=FALSE)
		n <- tabulate(z)
		if(all(n>=min.n)) return(list(group=z,breaks=breaks))
	}

#	Try equally spaced by quantiles
	z <- cut(x,breaks=breaks,labels=FALSE)
	n <- tabulate(z)
	if(all(n>=min.n)) return(list(group=z,breaks=breaks))

#	If all else fails, order by x
	o <- order(x)
	n <- floor(nx/intervals)
	nresid <- nx - intervals*n
	n <- rep_len(n,intervals)
	if(nresid>0) n[1:nresid] <- n[1:nresid]+1
	z <- rep(1:intervals,n)
	z[o] <- z
	return(list(group=z,breaks=breaks.eqn))

#	Function should never fail
	stop("Could not cut x into requested number of intervals with specified min.n in each group")
}
