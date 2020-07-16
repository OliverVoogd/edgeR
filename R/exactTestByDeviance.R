exactTestByDeviance <-  function(y1,y2,dispersion=0)
#	Test for differences in means between two groups of
#	negative binomial or Poisson random variables,
#	using exact enumeration conditional on total sum.

#	Rejection region is defined by large deviance statistics,
#	so this is a conditional likelihood ratio test.

#	R version by Davis McCarthy, Gordon Smyth, created 8 August 2011.
#	C++ version by Aaron Lun, created 26 July 2012.
#	Last modified 9 Dec 2013.
{
	y1 <- as.matrix(y1)
	y2 <- as.matrix(y2)
	ntags <- nrow(y1)
	if(ntags!=nrow(y2)) stop("Number of rows of y1 not equal to number of rows of y2")
	if(any(is.na(y1)) || any(is.na(y2))) stop("NAs not allowed")
	n1 <- ncol(y1)
	n2 <- ncol(y2)

	if(n1==n2) return(exactTestDoubleTail(y1=y1,y2=y2,dispersion=dispersion))

	sum1 <- as.integer(round(rowSums(y1)))
	sum2 <- as.integer(round(rowSums(y2)))
	if(all(dispersion==0)) return(binomTest(sum1,sum2,p=n1/(n1+n2)))
	if(any(dispersion==0)) stop("dispersion must be either all zero or all positive")
	if(length(dispersion)==1) dispersion <- rep(dispersion,ntags)

	pvals <- rep(1,ntags)
	if(ntags==0) return(pvals)

#	Eliminate all zero rows
	all.zeros <- sum1==0 & sum2==0
	if(any(all.zeros)) {
		pvals[!all.zeros] <- Recall(y1=y1[!all.zeros,,drop=FALSE],y2=y2[!all.zeros,,drop=FALSE],dispersion=dispersion[!all.zeros])
		return(pvals)
	}

#	Checking the dispersion type.
	dispersion<-as.double(dispersion)
	pvals<-.Call(.cxx_exact_test_by_deviance, sum1, sum2, n1, n2, dispersion)
	pmin(pvals, 1)
}
