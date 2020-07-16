# Calculate pseudo counts and pseudo library sizes.

equalizeLibSizes <- function(y, ...)
UseMethod("equalizeLibSizes")

equalizeLibSizes.DGEList <- function(y, dispersion=NULL, ...)
#	Yunshun Chen. Created 17 March 2016.
{
#	Check y
	y <- validDGEList(y)

#	Check dispersion
	if(is.null(dispersion)) dispersion <- getDispersion(y)
	
	lib.size <- y$samples$lib.size * y$samples$norm.factors

	out <- equalizeLibSizes(y=y$counts, group=y$samples$group, dispersion=dispersion, lib.size=lib.size)
	y$pseudo.counts <- out$pseudo.counts
	y$pseudo.lib.size <- out$pseudo.lib.size
	y
}

equalizeLibSizes.default <- function(y, group=NULL, dispersion=NULL, lib.size=NULL, ...)
#	Uses a quantile-to-quantile transformation so that new count counts are equivalent deviates on the equalized scale.
#	Davis McCarthy, Gordon Smyth.
#	Created July 2009. Last modified 17 March 2016.
{
#	Check y
	y <- as.matrix(y)
	ntags <- nrow(y)
	nlibs <- ncol(y)

#	Check group
	if(is.null(group)) group <- rep(1, nlibs)
	if(length(group)!=nlibs) stop("Incorrect length of group.")
	group <- dropEmptyLevels(group)

#	Check dispersion
	if(is.null(dispersion)) dispersion <- 0.05
	
#	Check lib.size
	if(is.null(lib.size)) lib.size <- colSums(y)
	if(length(lib.size)!=nlibs) stop("Incorrect length of lib.size.")
		
	common.lib.size <- exp(mean(log(lib.size)))
	levs.group <- unique(group)
	input.mean <- output.mean <- matrix(0, ntags, nlibs)
	for(i in 1:length(levs.group)) {
		j <- group==levs.group[i]
		beta <- mglmOneGroup(y[,j,drop=FALSE], dispersion=dispersion, offset=log(lib.size[j]))
		lambda <- exp(beta)
		input.mean[,j] <- matrix(lambda,ncol=1) %*% matrix(lib.size[j],nrow=1)
		output.mean[,j] <- matrix(lambda, ncol=1) %*% matrix(common.lib.size, nrow=1, ncol=sum(j))
	}
	pseudo <- q2qnbinom(y, input.mean=input.mean, output.mean=output.mean, dispersion=dispersion)
	pseudo[pseudo<0] <- 0
	list(pseudo.counts=pseudo, pseudo.lib.size=common.lib.size)
}
