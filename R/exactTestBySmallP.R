exactTestBySmallP <- function(y1,y2,dispersion=0) 
#	Test for differences in means between two groups of
#	negative binomial or Poisson random variables,
#	using exact enumeration conditional on total sum.

#	Rejection region is by method of small probability, i.e.,
#	all values with probability equal or less than that observed.

#	Mark Robinson, Davis McCarthy, Gordon Smyth.
#	Created 17 June 2009.  Last modified 9 Dec 2013.
{
	y1 <- as.matrix(y1)
	y2 <- as.matrix(y2)
	ntags <- nrow(y1)
	if(ntags!=nrow(y2)) stop("Number of rows of y1 not equal to number of rows of y2")
	if(any(is.na(y1)) || any(is.na(y2))) stop("NAs not allowed")
	n1 <- ncol(y1)
	n2 <- ncol(y2)

	if(n1==n2) return(exactTestDoubleTail(y1=y1,y2=y2,dispersion=dispersion))

	sum1 <- round(rowSums(y1))
	sum2 <- round(rowSums(y2))
	N <- sum1+sum2
	mu <- N/(n1+n2)
	if(all(dispersion==0)) return(binomTest(sum1,sum2,p=n1/(n1+n2)))
	if(any(dispersion==0)) stop("dispersion must be either all zero or all positive")
	if(length(dispersion)==1) dispersion <- rep(dispersion,ntags)
	r <- 1/dispersion
	all.zeros <- N==0

	pvals <- rep(1,ntags)
	if(ntags==0) return(pvals)
	if(any(all.zeros)) {
		pvals[!all.zeros] <- Recall(y1=y1[!all.zeros,,drop=FALSE],y2=y2[!all.zeros,,drop=FALSE],dispersion=dispersion[!all.zeros])
		return(pvals)
	}
	for (i in 1:ntags) {
		ind <- 0:N[i]
		p.top <- dnbinom(ind,size=n1*r[i],mu=n1*mu[i])*dnbinom(N[i]-ind,size=n2*r[i],mu=n2*mu[i])
		p.obs <- dnbinom(sum1[i],size=n1*r[i],mu=n1*mu[i]) * dnbinom(sum2[i],size=n2*r[i],mu=n2*mu[i])
		keep <-  p.top<=p.obs
		p.bot <- dnbinom(N[i],size=(n1+n2)*r[i],mu=(n1+n2)*mu[i])
		pvals[i] <- sum(p.top[keep]/p.bot)
	}
	min(pvals,1)
}
