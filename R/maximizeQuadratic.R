maximizeQuadratic <- function(y, x=1:ncol(y))
#	Maximize a function using quadratic approximation
#	Yunshun Chen and Gordon Smyth
#	Created 20 Oct 2011.  Last modified 7 Aug 2019.
{	
	if(is.vector(y)) 
		y <- matrix(y, nrow=1)
	if(!is.matrix(y)) stop("y must be either a vector or a matrix.")
	npts <- ncol(y)
	ntags <- nrow(y)
	i <- max.col(y)

	if(is.matrix(x)){
		if(!all(dim(x)==dim(y))) stop("The dimensions of x and y not matched.")
	} else if(is.vector(x)){
		if(length(x)!=npts) stop("The length of the vector x does not match the dimension of y.")
		x <- matrix(x, ntags, npts, byrow=TRUE)
	} else {
		stop("x must be either a vector or a matrix")
	}

	z <- rep_len(NA_real_, ntags)
	lower <- i==1
	upper <- i==npts
	z[lower] <- x[lower, 1]
	z[upper] <- x[upper, npts]

	keep <- !lower & !upper
	i <- i[keep]
	w <- y[keep, ]
	x <- x[keep, ]
	ntags <- length(i)
	
	r <- seq(from=0, by=npts, len=ntags)
	y1 <- t(w)[r+i-1]
	y2 <- t(w)[r+i]
	y3 <- t(w)[r+i+1]

	x1 <- t(x)[r+i-1]
	x2 <- t(x)[r+i]
	x3 <- t(x)[r+i+1]

	mid12 <- (x1+x2)/2
	mid23	<- (x2+x3)/2
	
	d21 <- (y2-y1)/(x2-x1)
	d32 <- (y3-y2)/(x3-x2)

	z[keep] <- (d21*mid23-d32*mid12)/(d21-d32)
	z
}

