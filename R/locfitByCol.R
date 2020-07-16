locfitByCol <- function(y, x=NULL, weights=1, span=0.5, degree=0)
#	Gordon Smyth
#	20 Aug 2012.  Last modified 15 June 2016.
{
	y <- as.matrix(y)
	ntags <- nrow(y)
	weights <- rep_len(weights,ntags)
	if(is.null(x)) x <- 1:ntags
	if(span*ntags<2 || ntags<=1) return(y)
	for (j in 1:ncol(y)) y[,j] <- fitted(locfit(y[,j]~x,weights=weights, alpha=span,deg=degree))
	y
}
