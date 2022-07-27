zscoreNBinom <- function(q, size, mu, method="midp")
#	Z-score equivalents for negative binomial deviates
#	Non-integer values for q are allowed
#	Gordon Smyth, Aaron Lun
#	Created 10 December 2011
#	Last modified 15 Feb 2019
{
#	Ensure arguments all same length
	n <- length(q)
	size <- rep_len(size,length.out=n)
	mu <- rep_len(mu,length.out=n)

#	Check method
	method <- match.arg(method, c("midp","random"))
	if(method=="randomized") return(.zscoreNBinomRandomized(q=q,size=size,mu=mu))

#	Output object
	z <- q

#	Point mass, allowing for non-integer values
	qr <- round(q)
	logd <- dnbinom(qr,size=size,mu=mu,log=TRUE)

#	Position of q relative to nearest integer
	delta <- q - qr

#	Treat q=0 as special case (not necessary, just for efficiency)
	Z <- (qr == 0)
	w <- delta[Z] + 0.5
	logp <- logd[Z] + log(w)
	z[Z] <- qnorm(logp,lower.tail=TRUE,log.p=TRUE)

#	Right and left tails
	R <- (q >= mu)
	L <- which(!R & !Z)
	R <- which(R & !Z)

#	Right tail scores
	if(length(R)) {
#		Tail prob w/o point mass
		logp <- pnbinom(qr[R]+0.5,size=size[R],mu=mu[R],lower.tail=FALSE,log.p=TRUE)

#		Weight for point mass
		w <- 0.5 - delta[R]

#		Tail prob with weighted point mass
		logp <- logsumexp(logp, logd[R] + log(w))

#		Z-score
		z[R] <- qnorm(logp,lower.tail=FALSE,log.p=TRUE)
	}

#	Left tail scores
	if(length(L)) {
#		Tail prob w/o point mass
		logp <- pnbinom(qr[L]-0.5,size=size[L],mu=mu[L],lower.tail=TRUE,log.p=TRUE)

#		Weight for point mass
		w <- delta[L] + 0.5

#		Tail prob with weighted point mass
		logp <- logsumexp(logp, logd[L] + log(w))

#		Z-score
		z[L] <- qnorm(logp,lower.tail=TRUE,log.p=TRUE)
	}

	z
}
