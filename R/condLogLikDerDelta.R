condLogLikDerDelta <- function(y,delta,der=1L)
# Derivatives of log-likelihood function wrt to delta
# r=1/dispersion and delta=1/(1+r)=dispersion/(1+dispersion)
# der is order of derivative required (0th deriv is the function)
# Written by Mark Robinson, edited by Davis McCarthy, February 2009
{
#	Vector interpreted as matrix of one row, i.e., one gene
	if (is.vector(y)) {
		y <- matrix(y,nrow=1)
	} else {
		y <- as.matrix(y)
	}
	if( !(length(delta)==1 | length(delta)==nrow(y)) ) stop("delta must be of length 1 or nrow(y)")

	r <- (1/delta)-1
	switch(der+1L,
		condLogLikDerSize(y,r,der=0L),
		condLogLikDerSize(y,r,der=1L)*(-delta^(-2)),
		condLogLikDerSize(y,r,der=1L)*2*(delta^(-3))+condLogLikDerSize(y,r,der=2)*(delta^(-4))
	)
}

