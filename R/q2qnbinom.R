q2qpois <- function (x, input.mean, output.mean)
#	Approximate quantile to quantile mapping between Poisson distributions
#	Original version, Gordon Smyth, 31 July 2009
{
	if(any(x<0)) stop("x must be non-negative")
	if(any(input.mean<0)) stop("input.mean must be non-negative")
	if(any(output.mean<0)) stop("output.mean must be non-negative")
	eps <- 1e-14
	zero <- input.mean<eps | output.mean<eps
	input.mean[zero] <- input.mean[zero]+0.25
	output.mean[zero] <- output.mean[zero]+0.25
	i <- (x >= input.mean)
	j <- !i
	p1 <- p2 <- q1 <- q2 <- x
	if(any(i)) {
		p1[i] <- pnorm(x[i], mean=input.mean[i], sd=sqrt(input.mean[i]), lower.tail=FALSE, log.p=TRUE)
		p2[i] <- pgamma(x[i], shape=(input.mean[i]), lower.tail=FALSE, log.p=TRUE)
		q1[i] <- qnorm(p1[i], mean=output.mean[i], sd=sqrt(output.mean[i]), lower.tail=FALSE, log.p=TRUE)
		q2[i] <- qgamma(p2[i], shape=(output.mean[i]), lower.tail=FALSE, log.p=TRUE)
	}
	if(any(j)) {
		p1[j] <- pnorm(x[j], mean=input.mean[j], sd=sqrt(input.mean[j]), lower.tail=TRUE, log.p=TRUE)
		p2[j] <- pgamma(x[j], shape=(input.mean[j]), lower.tail=TRUE, log.p=TRUE)
		q1[j] <- qnorm(p1[j], mean=output.mean[j], sd=sqrt(output.mean[j]), lower.tail=TRUE, log.p=TRUE)
		q2[j] <- qgamma(p2[j], shape=(output.mean[j]), lower.tail=TRUE, log.p=TRUE)
	}
	(q1+q2)/2        
}

q2qnbinom <- function(x, input.mean, output.mean, dispersion=0)
#	Approximate quantile to quantile mapping between negative-binomial distributions
#	with different means but same dispersion
#	Original version, Gordon Smyth, 31 July 2009
{
	if(any(x<0)) stop("x must be non-negative")
	if(any(input.mean<0)) stop("input.mean must be non-negative")
	if(any(output.mean<0)) stop("output.mean must be non-negative")
	if(any(dispersion<0)) stop("dispersion must be non-negative")
	eps <- 1e-14
	zero <- input.mean<eps | output.mean<eps
	input.mean[zero] <- input.mean[zero]+0.25
	output.mean[zero] <- output.mean[zero]+0.25
	ri <- 1+dispersion*input.mean
	vi <- input.mean*ri
	ro <- 1+dispersion*output.mean
	vo <- output.mean*ro
	i <- (x >= input.mean)
	j <- !i
	p1 <- p2 <- q1 <- q2 <- x
	if(any(i)) {
		p1[i] <- pnorm(x[i], mean=input.mean[i], sd=sqrt(vi[i]), lower.tail=FALSE, log.p=TRUE)
		p2[i] <- pgamma(x[i], shape=input.mean[i]/ri[i], scale=ri[i], lower.tail=FALSE, log.p=TRUE)
		q1[i] <- qnorm(p1[i], mean=output.mean[i], sd=sqrt(vo[i]), lower.tail=FALSE, log.p=TRUE)
		q2[i] <- qgamma(p2[i], shape=output.mean[i]/ro[i], scale=ro[i], lower.tail=FALSE, log.p=TRUE)
	}
	if(any(j)) {
		p1[j] <- pnorm(x[j], mean=input.mean[j], sd=sqrt(vi[j]), lower.tail=TRUE, log.p=TRUE)
		p2[j] <- pgamma(x[j], shape=input.mean[j]/ri[j], scale=ri[j], lower.tail=TRUE, log.p=TRUE)
		q1[j] <- qnorm(p1[j], mean=output.mean[j], sd=sqrt(vo[j]), lower.tail=TRUE, log.p=TRUE)
		q2[j] <- qgamma(p2[j], shape=output.mean[j]/ro[j], scale=ro[j], lower.tail=TRUE, log.p=TRUE)
	}
	(q1+q2)/2
}
