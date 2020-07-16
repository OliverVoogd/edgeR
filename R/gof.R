gof <- function( glmfit, pcutoff=0.1, adjust="holm", plot=FALSE, main="qq-plot of residual deviances", ...)
#	Compare residual deviances to chisquare distribution to identify dispersion outlier genes
#	Davis McCarthy, Gordon Smyth
#	23 Mar 2011. Last modified 2 Jun 2020.
{
	stopifnot( is(glmfit, "DGEGLM") )
	gof.stats <- glmfit$deviance
	gof.pvals <- pchisq(gof.stats, df=glmfit$df.residual, lower.tail=FALSE, log.p=FALSE)
	outlier <- p.adjust(gof.pvals, method=adjust) < pcutoff

	if(plot) {
		n <- length(gof.stats)
		z <- zscoreGamma(gof.stats,shape=glmfit$df.residual/2,scale=2)
		col <- rep_len("black",n)
		col[outlier] <- "blue"
		pch <- rep_len(1,n)
		pch[outlier] <- 16
		qqnorm(z,col=col,pch=pch,main=main,...)
		abline(0,1)
	}

	invisible(list(gof.statistics=gof.stats, gof.pvalues=gof.pvals, outlier=outlier, df=glmfit$df.residual[1]))
}
