goodTuring <- function(x, conf=1.96)
#	Simple Good-Turing algorithm for frequency estimation
#	as described by Gale and Sampson (1995)

#  Sampson's C code translated to C++ and optimized by Aaron Lun
#	Has been tested against Sampson's C code from
#	http://www.grsampson.net/RGoodTur.html
#	and gives identical results.

#	Gordon Smyth and Aaron Lun
#	9 Nov 2010.  Last modified 7 Sep 2012.
{
#	Raw frequencies
	x <- as.integer(x)
	if(max(x) < length(x)) {
		n <- tabulate(x+1L)
		n0 <- n[1]
		n <- n[-1]
		max.x <- length(n)
		r <- seq.int(from=1L,to=max.x)
		r <- r[n>0]
		n <- n[n>0]
	} else {
		r <- sort(unique(x))
		z <- match(x,r)
		n <- tabulate(z)
		if(r[1]==0) {
			n0 <- n[1]
			n <- n[-1]
			r <- r[-1]
		} else {
			n0 <- 0
		}
	}

#	SGT algorithm (no type checking, as that's enforced above)
	out <- .Call(.cxx_simple_good_turing, r, n, conf)
	names(out) <- c("P0","proportion")

	out$count <- r
	out$n <- n
	out$n0 <- n0
	out
}

goodTuringPlot <- function(x)
#	Simple Good-Turing algorithm for frequency estimation
#	as described by Gale and Sampson (1995)

#	Has been tested against Sampson's C code from
#	http://www.grsampson.net/RGoodTur.html
#	and gives identical results.

#	Gordon Smyth
#	9 Nov 2010.  Last modified 21 Nov 2010.
{
#	Raw frequencies
	n <- tabulate(x+1L)
	n0 <- n[1]
	n <- n[-1]
	max.x <- length(n)
	r <- seq.int(from=1L,to=max.x)

#	Fit a linear trend to log-frequencies
	n.pos <- n[n>0]
	r.pos <- r[n>0]
	l <- length(n.pos)
	q <- diff(c(0L,r.pos,2L*r.pos[l]-r.pos[l-1]),lag=2)/2
	z <- n.pos/q
	design <- cbind(1,log(r.pos))
	fit <- lm.fit(x=design,y=log(z))
	plot(log(r.pos),log(z),xlab="log count",ylab="log frequency")
	abline(fit)
	invisible(list(r=r.pos,n=n.pos))
}

goodTuringProportions <- function(counts)
#	Transform counts to approximately normal expression values
#	Gordon Smyth
#	15 Dec 2010.  Last modified 5 Jan 2011.
{
	z <- counts <- as.matrix(counts)
	nlibs <- ncol(counts)
	for (i in 1:nlibs) {
		g <- goodTuring(counts[,i])
		p0 <- g$P0/g$n0
		zero <- z[,i]==0
		z[zero,i] <- p0
		m <- match(z[!zero,i],g$count)
		z[!zero,i] <- g$proportion[m]
	}
	z
}
