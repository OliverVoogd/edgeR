dropEmptyLevels <- function(x)
#	Drop levels of a factor that don't occur
#	Gordon Smyth
#	Created 25 March 2012.  Last modified 6 March 2015.
{
	if(is.factor(x)) {
		i <- which(tabulate(as.integer(x))>0L)
		if(length(i) < nlevels(x)) x <- factor(x, levels=levels(x)[i])
		return(x)
	} else {
		return(factor(x))
	}
}
