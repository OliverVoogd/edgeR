weightedCondLogLikDerDelta <- function(y, delta, tag, prior.n=10, ntags=nrow(y[[1]]), der=0)
# Weighted conditional log-likelihood for a tag - necessary to estimate tagwise dispersions
# Created 10 Sep 2009. Last modified 2 Jun 2020.
{
	l0 <- rep_len(0,ntags)
	onev <- rep_len(1,ntags)
	for(i in seq_len(length(y))) {
		l0 <- condLogLikDerDelta(y[[i]],delta,der=der)+l0
	}
	m0 <- sum(l0)
	l0a <- l0 + (prior.n/ntags)*m0
	l0a[tag]
}
