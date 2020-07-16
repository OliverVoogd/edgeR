calcNormOffsetsforChIP <- function(input,response,dispersion=0.01,niter=6,loss="p",plot=FALSE,verbose=FALSE,...)
#	Normalize ChIP-Seq counts to input and form offset matrix
#	Gordon Smyth  
#	14 Dec 2011.  Last modified 14 May 2012.
{
	input <- as.matrix(input)
	y <- as.matrix(response)

	if(nrow(input) != nrow(y)) stop("nrows of input and response disagree")
	if(ncol(input)==1 && ncol(y)>1) input <- matrix(input,nrow(input),ncol(response))
	if(ncol(input) != ncol(y)) stop("ncols of input and response disagree")

	offset <- y
	for (j in 1:ncol(y)) {
		out <- normalizeChIPtoInput(input[,j],y[,j],dispersion=dispersion,niter=niter,loss=loss,plot=plot,verbose=verbose,main=colnames(y)[j],...)
		offset[,j] <- log(out$scaling.factor * input[,j])
	}

	if(is(response,"DGEList")) {
		response$offset <- offset
		return(response)
	} else {
		return(offset)
	}
}
