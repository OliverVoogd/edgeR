movingAverageByCol <- function(x,width=5,full.length=TRUE)
#	Moving average smoother for columns of a matrix
#	Gordon Smyth
#	17 Feb 2011
{
	x <- as.matrix(x)
	width <- as.integer(width)
	if(width<=1) return(x)
	n <- nrow(x)
	m <- ncol(x)
	if(width>n) {
		width <- n
		warning("reducing moving average width to nrow(x)")
	}
	if(full.length) {
		half1 <- ceiling(width/2)
		half2 <- floor(width/2)
		x <- rbind(matrix(0,half1,m),x,matrix(0,half2,m))
	} else {
		if(width==n) return(matrix(colMeans(x),1L,m))
		x <- rbind(matrix(0,1,m),x)
	}
	n2 <- nrow(x)
	x <- apply(x,2,cumsum)
	x <- x[(width+1):n2,,drop=FALSE]-x[1:(n2-width),,drop=FALSE]
	n3 <- nrow(x)
	w <- rep(width,n3)
	if(full.length) {
		if(half1>1) w[1:(half1-1)] <- width-(half1-1):1
		w[(n3-half2+1):n3] <- width-(1:half2)
	}
	x/w
}

