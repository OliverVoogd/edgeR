#  FIT GENERALIZED LINEAR MODELS

filterByExpr <- function(y, ...)
UseMethod("filterByExpr")

filterByExpr.DGEList <- function(y, design=NULL, group=NULL, lib.size=NULL, ...)
{
#	Order of precedence:
#	1. group as argument
#	2. design as argument
#	3. y$design
#	4. y$samples$group
	if(is.null(design) && is.null(group)) {
		design <- y$design
		if(is.null(design)) {
			group <- y$samples$group
			if(length(levels(group))==1L) warning("All samples appear to belong to the same group.")
		}
	}
	if(is.null(lib.size)) lib.size <- y$samples$lib.size * y$samples$norm.factors
	filterByExpr.default(y$counts, design=design, group=group, lib.size=lib.size, ...)
}

filterByExpr.SummarizedExperiment <- function(y, design=NULL, group=NULL, lib.size=NULL, ...)
#	Created 19 March 2020. Last revised 19 March 2020.
{
	y <- SE2DGEList(y)
	filterByExpr.DGEList(y, design=design, group=group, lib.size=lib.size, ...)
}

filterByExpr.default <- function(y, design=NULL, group=NULL, lib.size=NULL, min.count=10, min.total.count=15, large.n=10, min.prop=0.7, ...)
#	Filter low expressed genes given count matrix
#	Computes TRUE/FALSE index vector indicating which rows to keep
#	Gordon Smyth
#	Created 13 Nov 2017. Last revised 26 Jan 2020.
{
	y <- as.matrix(y)
	if(mode(y) != "numeric") stop("y is not a numeric matrix")
	if(is.null(lib.size)) lib.size <- colSums(y)

#	Minimum effect sample sample size for any of the coefficients
	if(is.null(group)) {
		if(is.null(design)) {
			message("No group or design set. Assuming all samples belong to one group.")
			MinSampleSize <- ncol(y)
		} else {
			h <- hat(design)
			MinSampleSize <- 1/max(h)
		}
	} else {
		group <- as.factor(group)
		n <- tabulate(group)
		MinSampleSize <- min(n[n > 0L])
	}
	if(MinSampleSize > large.n) MinSampleSize <- large.n + (MinSampleSize-large.n)*min.prop

#	CPM cutoff
	MedianLibSize <- median(lib.size)
	CPM.Cutoff <- min.count/MedianLibSize*1e6
	CPM <- cpm(y,lib.size=lib.size)
	tol <- 1e-14
	keep.CPM <- rowSums(CPM >= CPM.Cutoff) >= (MinSampleSize - tol)

#	Total count cutoff
	keep.TotalCount <- (rowSums(y) >= min.total.count - tol)

	keep.CPM & keep.TotalCount
}
