commonCondLogLikDerDelta <- function(y, delta, der=0) 
# Calculates the common conditional log-likelihood (i.e. summed over all tags) - necessary so that optimize can be applied in estimateCommonDisp
# Davis McCarthy, July 2009
{
	l0 <- 0
	for(i in 1:length(y)) {
		l0 <- condLogLikDerDelta(y[[i]],delta,der=der)+l0
	}
	sum(l0)
}
