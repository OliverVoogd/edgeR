plotMeanVar2 <- function(y, design=NULL, dispersion=0, offset=0, nbins=100, make.plot=TRUE, xlab="Mean", ylab="Ave. binned standardized residual", ... )
# Fit a Poisson (or NB) glm to each row
# Compute squared residuals
# Bin by size of fitted value and plot average squared residual vs average fitted value
# Davis McCarthy and Gordon Smyth
# Created 9 November 2010.
# Renamed from binStdResidPois() to dglmStdResid() and moved from meanVar.R to dglmStdResid.R on 25 Nov 2010.
# Renamed from dglmStdResid() to meanVar2() on 7 Aug 2019.
# Last modified 7 Aug 2019.
{
#	Fit Poisson or NB GLM
	ngenes <- nrow(y)
	nlibs <- ncol(y)
	if(is.null(design)) design <- matrix(1,nlibs,1)
	Mu <- glmFit(y, design=design, dispersion=0, offset=offset)$fitted.values

#	Extract fitted values and squared residuals
#	Use approximate leverages
	ColMeanMu <- colMeans(Mu)
	V <- ColMeanMu + dispersion*ColMeanMu^2
	h1 <- 1-hat(design/sqrt(V),intercept=FALSE)
	j <- (h1 < 1e-14)
	if(any(j)) {
		y <- y[,!j]
		Mu <- Mu[,!j]
		h1 <- h1[!j]
	}
	ResidSqr <- (y - Mu)^2 / h1

#	Bin by Mu and compute average Mu and ResidSqr for each bin
	bins <- cutWithMinN(Mu, intervals=nbins, min.n=3)
	X <- cbind(1,as.vector(Mu),as.vector(ResidSqr))
	BinMeans <- rowsum(X,bins$group,reorder=FALSE)
	AveMu <- BinMeans[,2] / BinMeans[,1]
	AveResidSqr <- BinMeans[,3] / BinMeans[,1]

#	Plot
	if(make.plot) plot(AveMu, AveResidSqr, log="xy", xlab=xlab, ylab=ylab, ...)

	invisible(list(mean=AveMu,var=AveResidSqr))
}
