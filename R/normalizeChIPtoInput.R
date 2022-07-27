normalizeChIPtoInput <- function(input,response,dispersion=0.01,niter=6,loss="p",plot=FALSE,verbose=FALSE,...)
#	Normalize ChIP-Seq counts to input
#	and test for enrichment
#	Gordon Smyth  
#	2 Dec 2011.  Last modified 1 Jun 2020.
{
	if(length(input)!=length(response)) stop("input and response must be same length")
	if(any(input<0) || any(response<0)) stop("negative values not allowed")
	if(any(dispersion<=0)) stop("dispersion must be positive")

#	Remove zero inputs from main calculation
	zero <- input<=0 & response<=0
	if(any(zero)) {
		p.value <- rep.int(1,length(zero))
		out <- Recall(input[!zero],response[!zero],dispersion=dispersion,niter=niter,loss=loss,plot=plot,verbose=verbose,...)
		p.value[!zero] <- out$p.value
		out$p.value <- p.value
		return(out)
	}

  	n <- length(response)

#	Handle special cases
  	if(n==0) return(list(p=numeric(0),scaling.factor=NA,prop.enriched=NA))
	if(all(input==0)) return(list(p=rep(0,1),scaling.factor=0,prop.enriched=1))
  	if(n==1) return(list(p=1,scaling.factor=input/response,prop.enriched=0))

#	Reset zero inputs to minimum positive value
	input[input==0] <- min(input[input>0])

#	From here, all values of input are positive

#	Objective function for optimizing scaling of response relative to input
	loss <- match.arg(loss,c("p","z"))
	f <- switch(loss,
		p = function(scaling.factor,input,response,prop.enriched) {
			p <- pnbinom(response,mu=scaling.factor*input,size=1/dispersion)
			d <- dnbinom(response,mu=scaling.factor*input,size=1/dispersion)
			pmid <- p-d/2
			n <- length(response)
			n.not.enriched <- round(length(response) * (1-prop.enriched))
			n.not.enriched <- max(n.not.enriched,1)
			p.sorted <- sort(pmid,partial=n.not.enriched)
			out <- abs(mean(p.sorted[1:n.not.enriched])-0.5)
		},
		z = function(scaling.factor,input,response,prop.enriched) {
			z <- zscoreNBinom(response,mu=scaling.factor*input,size=1/dispersion)
			n <- length(response)
			n.not.enriched <- round(length(response) * (1-prop.enriched))
			n.not.enriched <- max(n.not.enriched,1)
			z.sorted <- sort(z,partial=n.not.enriched)
			out <- mean(abs(z.sorted[1:n.not.enriched]))
		}
	)

#	Starting value for proportion of enriched marks
	prop.enriched <- 0.5
	scaling.factor.interval <- quantile(response/input,probs=c(0.1,0.8))
	
	if(diff(scaling.factor.interval)==0) {
		scaling.factor <- scaling.factor.interval[1]
		p <- pnbinom(response,mu=scaling.factor*input,size=1/dispersion,lower.tail=FALSE)
		d <- dnbinom(response,mu=scaling.factor*input,size=1/dispersion)
		pmid <- p-d/2
		enriched <- p.adjust(pmid,method="holm")<0.5
		prop.enriched <- sum(enriched)/n
		if(verbose) cat("prop.enriched:",prop.enriched,"scaling.factor:",scaling.factor,"\n")
	} else {
#		Iterate over prop.enriched and scaling.factor
		for (iter in 1:niter) {
			scaling.factor <- optimize(f,interval=scaling.factor.interval,input=input,response=response,prop.enriched=prop.enriched)$minimum
			p <- pnbinom(response,mu=scaling.factor*input,size=1/dispersion,lower.tail=FALSE)
			d <- dnbinom(response,mu=scaling.factor*input,size=1/dispersion)
			pmid <- p-d/2
			enriched <- p.adjust(pmid,method="holm")<0.5
			prop.enriched <- sum(enriched)/n
			if(verbose) cat("prop.enriched:",prop.enriched,"scaling.factor:",scaling.factor,"\n")
		}
	}

	if(plot) {
		x <- log2(input)
		y <- log2(response)
		mu <- log2(input*scaling.factor)
		o <- order(x)
		plot(x,y,xlab="log2(Input)",ylab="log2(Response)",type="n",...)
		points(x[!enriched],y[!enriched],pch=16,cex=0.1)
		points(x[enriched],y[enriched],pch=16,cex=0.2,col="red")
		lines(x[o],mu[o],col="blue",lwd=2)
		legend("topleft",pch=16,col=c("red","blue"),legend=c("Enriched","Line of normalization"))
	}

	return(list(p.value=p,pmid.value=pmid,scaling.factor=scaling.factor,prop.enriched=prop.enriched))
}
