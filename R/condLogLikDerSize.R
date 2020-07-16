condLogLikDerSize <- function(y, r, der=1L)
#	Derivatives of the conditional log-likelihood function (given the row sum)
#	with respect to r=1/dispersion
#	for a single group of replicate libraries, all of the same total size
{
#	Vector interpreted as matrix of one row, i.e., one gene
	if (is.vector(y)) {
		y <- matrix(y,nrow=1)
	} else {
		y <- as.matrix(y)
	}

	n <- ncol(y)
	m <- rowMeans(y)

	switch(der+1L,
		rowSums(lgamma(y+r)) + lgamma(n*r) - lgamma(n*(m+r)) - n*lgamma(r),
		rowSums(digamma(y+r)) + n*digamma(n*r) - n*digamma(n*(m+r)) - n*digamma(r),
		rowSums(trigamma(y+r)) + n^2*trigamma(n*r) - n^2*trigamma(n*(m+r)) - n*trigamma(r)
	)
}
