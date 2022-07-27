expandAsMatrix <- function(x,dim=NULL,byrow=TRUE)
#	Convert scalar, row or column vector, or matrix, to be a matrix
#	Gordon Smyth, Yunshun Chen, Aaron Lun
#	created 26 Jan 2011.  
#   last modified 29 Jun 2016.
{
#	Check dim argument
	if(is.null(dim)) return(as.matrix(x))
	if(length(dim)<2) stop("dim must be numeric vector of length 2")
	dim <- round(dim[1:2])
	if(any(dim<0)) stop("negative dimensions not allowed")

#	x is a vector
	dx <- dim(x)
	if(is.null(dx)) {
		lx <- length(x)
		if(lx==1) return(matrix(x,dim[1],dim[2]))
		if(lx==dim[1] & lx==dim[2]) return(matrix(x,dim[1],dim[2],byrow=byrow))
		if(lx==dim[2]) return(matrix(x,dim[1],dim[2],byrow=TRUE))
		if(lx==dim[1]) return(matrix(x,dim[1],dim[2],byrow=FALSE))
		stop("x of unexpected length")
	}

#	x is a matrix or data.frame
	if(length(dx)<2) stop("x has less than 2 dimensions")
	if(length(dx)>2) stop("x has more than 2 dimensions")
	if(all(dx==dim)) return(as.matrix(x))

#	x is a CompressedMatrix
	if(is(x, "CompressedMatrix")) {
        return(Recall(as.matrix(x), dim=dim, byrow=byrow))
    }
	stop("x is matrix of wrong size")
}

