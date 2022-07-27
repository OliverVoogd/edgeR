nearestReftoX <- function(x, reference, ...)
#	Find nearest element of reference for each element of x
#	reference should be sorted in increasing order
#	Gordon Smyth
#	3 Jan 2018. Last modified 5 Jan 2018.
{
	nref <- length(reference)
	midpt <- (reference[-nref]+reference[-1L])/2
	findInterval(x,midpt,...)+1L
}
