exactTestDoubleTail <- function(y1,y2,dispersion=0,big.count=900)
#	Test for differences in means between two groups of
#	negative binomial or Poisson random variables,
#	using exact enumeration conditional on total sum.

#	Smaller tail probability is doubled to get p-value.
#	QUESTION: should we use sign(logFC) to choose which tail to evaluate
#	instead of trying to find smaller of tail probabilities?

#	Gordon Smyth
#	28 Sep 2019.  Last modified 2 Jun 2020.
{
#	Convert matrices to vectors
	ntags <- NROW(y1)
	n1 <- NCOL(y1)
	n2 <- NCOL(y2)
	if(n1>1) s1 <- round(rowSums(y1)) else s1 <- round(y1)
	if(n2>1) s2 <- round(rowSums(y2)) else s2 <- round(y2)
	if(length(dispersion)==1) dispersion <- rep_len(dispersion,ntags)

#	Null fitted values
	s <- s1+s2
	mu <- s/(n1+n2)
	mu1 <- n1*mu
	mu2 <- n2*mu

	pvals <- rep_len(1,ntags)
	names(pvals) <- names(y1)

#	Poisson case
	pois <- dispersion<=0
#	BINOMTEST DOESN'T USE EQUAL TAILED REJECTION REGION
	if(any(pois)) pvals[pois] <- binomTest(s1[pois],s2[pois],p=n1/(n1+n2))

#	Use beta approximation for large counts
	big <- s1>big.count & s2>big.count
	if(any(big)) {
		y1 <- as.matrix(y1)
		y2 <- as.matrix(y2)
		pvals[big] <- exactTestBetaApprox(y1[big,,drop=FALSE],y2[big,,drop=FALSE],dispersion[big])
	}

	p.bot <- size1 <- size2 <- rep_len(0,ntags)
	left <- s1<mu1 & !pois & !big
	if(any(left)) {
		p.bot[left] <- dnbinom(s[left],size=(n1+n2)/dispersion[left],mu=s[left])
		size1[left] <- n1/dispersion[left]
		size2[left] <- n2/dispersion[left]
		for (g in which(left)) {
			x <- 0:s1[g]
			p.top <- dnbinom(x,size=size1[g],mu=mu1[g]) * dnbinom(s[g]-x,size=size2[g],mu=mu2[g])
			pvals[g] <- 2*sum(p.top)
		}
		pvals[left] <- pvals[left]/p.bot[left]
	}
	right <- s1>mu1 & !pois & !big
	if(any(right)) {
		p.bot[right] <- dnbinom(s[right],size=(n1+n2)/dispersion[right],mu=s[right])
		size1[right] <- n1/dispersion[right]
		size2[right] <- n2/dispersion[right]
		for (g in which(right)) {
			x <- s1[g]:s[g]
			p.top <- dnbinom(x,size=size1[g],mu=mu1[g]) * dnbinom(s[g]-x,size=size2[g],mu=mu2[g])
			pvals[g] <- 2*sum(p.top)
		}
		pvals[right] <- pvals[right]/p.bot[right]
	}
	pmin(pvals,1)
}
