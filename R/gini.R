gini <- function(x)
#	Gini diversity index for columns of a numeric matrix
#	Gordon Smyth
#	Created 5 Feb 2016. Last revised 17 April 2017.
{
	x <- as.matrix(x)
	d <- dim(x)
	for (j in 1:d[2]) x[,j] <- sort.int(x[,j],na.last=TRUE)
	i <- 1:d[1]
	m <- 0.75*d[1]
	S1 <- colSums((i-m)*x , na.rm=TRUE)
	S2 <- colSums(x, na.rm=TRUE)
	(2*(S1/S2+m)-d[1]-1)/d[1]
}

