exactTestBetaApprox <- function(y1,y2,dispersion=0)
#	Test for differences in means between two negative binomial
#	or Poisson random variables, or between two groups of variables,
#	using a beta distribution approximation.
#	Test is naturally conditional on total sum.
#	Left and right rejection regions have equal probability.

#	Gordon Smyth
#	28 Sep 2019.  Last modified 28 Sep 2011.
{
#	Convert matrices to vectors
	ntags <- NROW(y1)
	n1 <- NCOL(y1)
	n2 <- NCOL(y2)
	if(n1>1) y1 <- rowSums(y1)
	if(n2>1) y2 <- rowSums(y2)
	if(length(dispersion)==1) dispersion <- rep(dispersion,ntags)

#	Null fitted values
	y <- y1+y2
	mu <- y/(n1+n2)

#	Compute p-values
	pvals <- rep(1,ntags)
	all.zero <- y<=0
	alpha1 <- n1*mu/(1+dispersion*mu)
	alpha2 <- n2/n1*alpha1
	med <- rep(0,ntags)
	med[!all.zero] <- qbeta(0.5,alpha1[!all.zero],alpha2[!all.zero])
	left <- (y1+0.5)/y<med & !all.zero
	if(any(left)) {
		pvals[left] <- 2*pbeta((y1[left]+0.5)/y[left],alpha1[left],alpha2[left])
	}
	right <- (y1-0.5)/y>med & !all.zero
	if(any(right)) {
		pvals[right] <- 2*pbeta((y1[right]-0.5)/y[right],alpha1[right],alpha2[right],lower.tail=FALSE)
	}
	names(pvals) <- names(y1)
	pvals
}
