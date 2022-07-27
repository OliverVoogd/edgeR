systematicSubset <- function(n,order.by)
#	Take a systematic subset of indices,
#	stratified by a ranking variable
#	Gordon Smyth
#	28 Jan 2011
{
	ntotal <- length(order.by)
	sampling.ratio <- floor(ntotal/n)
	if(sampling.ratio <= 1) return(1:ntotal)
	i1 <- floor(sampling.ratio/2)+1
	i <- seq.int(from=i1,to=ntotal,by=sampling.ratio)
	o <- order(order.by)
	o[i]
}
