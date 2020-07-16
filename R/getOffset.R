getOffset <- function(y)
#	Extract offset vector or matrix from data object and optional arguments.
#	By default, offset is constructed from the lib.size and norm.factors
#	but offset supplied explicitly takes precedence

#	Gordon Smyth
#	26 Jan 2011. Last modified 11 Jan 2012.
{
	offset <- y$offset
	lib.size <- y$samples$lib.size
	norm.factors <- y$samples$norm.factors
	
	if(!is.null(offset)) {
		return(offset)
	} else {		
		if(!is.null(norm.factors)) lib.size <- lib.size*norm.factors
		return(log(lib.size))
	}
}
