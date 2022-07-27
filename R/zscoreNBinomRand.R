.zscoreNBinomRandomized <- function(q, size, mu)
#	Negative binomial randomized quantile residuals
#	Assumes integer values for q
#	Gordon Smyth
#	Created 15 Feb 2019. Last modified 15 Feb 2019.
{
	n <- length(q)
	q <- round(q)

#	Output object
	z <- q

#	Right and left tails
	R <- (q >= mu)
	L <- which(!R)
	R <- which(R)

#	Right tail scores
	if(length(R)) {
#		Tail prob with point mass
		p.with <- pnbinom(q[R]-0.5,size=size[R],mu=mu[R],lower.tail=FALSE,log.p=FALSE)

#		Tail prob without point mass
		p.without <- pnbinom(q[R]+0.5,size=size[R],mu=mu[R],lower.tail=FALSE,log.p=FALSE)

#		Randomize p
		p <- runif(length(R),min=p.without,max=p.with)

#		Z-score
		z[R] <- qnorm(p,lower.tail=FALSE,log.p=FALSE)

#		Check for floating underflow of tail probability
		if(is.infinite(max(z[R]))) {
			i <- (is.infinite(z[R]) & is.finite(q[R]))
			if(any(i)) {
				j <- R[i]
				z[j] <- zscoreNBinom(q[j],size=size[j],mu=mu[j])
			}
		}
	}

#	Left tail scores
	if(length(L)) {
#		Tail prob with point mass
		p.with <- pnbinom(q[L]+0.5,size=size[R],mu=mu[R],lower.tail=TRUE,log.p=FALSE)

#		Tail prob without point mass
		p.without <- pnbinom(q[L]-0.5,size=size[R],mu=mu[R],lower.tail=TRUE,log.p=FALSE)

#		Randomize p
		p <- runif(length(L),min=p.without,max=p.with)

#		Z-score
		z[L] <- qnorm(p,lower.tail=TRUE,log.p=FALSE)
	}

	z
}
