## binomTest.R

binomTest <- function(y1, y2, n1=sum(y1), n2=sum(y2), p=n1/(n1+n2))
#	Multiple exact binomial tests.
#	Intended for comparing DGE libraries
#	in the absence of biological variation.
#	Rejection region is all values with lower prob than that of
#	value observed, same as for binom.test() in stats package.

#	Based on function sage.test() in the statmod package.
#	Gordon Smyth
#	In statmod package 15 Nov 2003.
#	In edgeR package 11 Feb 2011.
#	Last modified 9 Jun 2020.
{
	if(length(y1) != length(y2)) stop("y1 and y2 must have same length")
	if(any(is.na(y1)) || any(is.na(y2))) stop("missing values not allowed")
	y1 <- round(y1)
	y2 <- round(y2)
	if(any(y1<0) || any(y2<0)) stop("y1 and y2 must be non-negative")
	if(p<=0 || p>=1) stop("p must be between 0 and 1")
	size <- y1+y2
	p.value <- rep_len(1,length(y1))
	if(p==0.5) {
		i <- (size>0)
		if(any(i)) {
			y1 <- pmin(y1[i],y2[i])
			size <- size[i]
			p.value[i] <- pmin(2*pbinom(y1,size=size,prob=0.5),1)
		}
		return(p.value)
	}
	if(any(big <- size>10000)) {
		ibig <- which(big)
		for (i in ibig) p.value[i] <- chisq.test(matrix(c(y1[i],y2[i],n1-y1[i],n2-y2[i]),2,2))$p.value
	}
	size0 <- size[size>0 & !big]
	if(length(size0)) for (isize in unique(size0)) {
		i <- (size==isize)
		d <- dbinom(0:isize,prob=p,size=isize)
		o <- order(d)
		cumsump <- cumsum(d[o])[order(o)]
		p.value[i] <- cumsump[y1[i]+1]
	}
	p.value
}
