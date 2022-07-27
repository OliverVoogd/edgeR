roast.DGEList <- function(y, index=NULL, design=NULL, contrast=ncol(design), geneid=NULL, set.statistic="mean", gene.weights=NULL, nrot=1999, ...)
#	Rotation gene set testing for RNA-Seq data
#	Yunshun Chen, Gordon Smyth
#	Created 19 Dec 2012. Last revised on 28 July 2019.
{
	if(is.null(design)) {
		design <- y$design
		if(is.null(design)) design <- model.matrix(~y$samples$group)
	}
	if(length(geneid)==1L) geneid <- y$genes[[geneid]]
	y <- .zscoreDGE(y=y, design=design, contrast=contrast)
	roast(y=y, index=index, design=design, contrast=contrast, geneid=geneid, set.statistic=set.statistic, gene.weights=gene.weights, var.prior=1, df.prior=Inf, nrot=nrot, approx.zscore=TRUE)
}

mroast.DGEList <- function(y, index=NULL, design=NULL, contrast=ncol(design), geneid=NULL, set.statistic="mean", gene.weights=NULL, nrot=1999, adjust.method="BH", midp=TRUE, sort="directional", ...)
#	Rotation gene set testing for RNA-Seq data with multiple sets
#	Yunshun Chen, Gordon Smyth
#	Created 8 Jan 2013.  Last revised 28 July 2019.
{
	if(is.null(design)) {
		design <- y$design
		if(is.null(design)) design <- model.matrix(~y$samples$group)
	}
	if(length(geneid)==1L) geneid <- y$genes[[geneid]]
	y <- .zscoreDGE(y=y, design=design, contrast=contrast)
	mroast(y=y, index=index, design=design, contrast=contrast, geneid=geneid, set.statistic=set.statistic, gene.weights=gene.weights, var.prior=1, df.prior=Inf, nrot=nrot, approx.zscore=TRUE, adjust.method=adjust.method, midp=midp, sort=sort)
}

fry.DGEList <- function(y, index=NULL, design=NULL, contrast=ncol(design), geneid=NULL, sort="directional", ...)
#	Rotation gene set testing for RNA-Seq data with multiple sets
#	Yunshun Chen, Gordon Smyth
#	Created 1 Dec 2015.  Last revised 3 June 2018.
{
	if(is.null(design)) {
		design <- y$design
		if(is.null(design)) design <- model.matrix(~y$samples$group)
	}
	if(length(geneid)==1L) geneid <- y$genes[[geneid]]
	y <- .zscoreDGE(y=y, design=design, contrast=contrast)
	fry(y=y, index=index, design=design, contrast=contrast, standardize="none", geneid=geneid, sort=sort)
}

camera.DGEList <- function(y, index, design=NULL, contrast=ncol(design), weights=NULL, use.ranks=FALSE, allow.neg.cor=FALSE, inter.gene.cor=0.01, sort=TRUE, ...)
#	Rotation gene set testing for RNA-Seq data accounting for inter-gene correlation
#	Yunshun Chen, Gordon Smyth
#	Created 07 Jan 2013. Last modified 29 Dec 2017.
{
	if(is.null(design)) {
		design <- y$design
		if(is.null(design)) design <- model.matrix(~y$samples$group)
	}
	y <- .zscoreDGE(y=y, design=design, contrast=contrast)
	camera(y=y, index=index, design=design, contrast=contrast, weights=weights, use.ranks=use.ranks, allow.neg.cor=allow.neg.cor, inter.gene.cor=inter.gene.cor, trend.var=FALSE, sort=sort)
}

romer.DGEList <- function(y, index, design=NULL, contrast=ncol(design), ...)
#	rotation mean-rank version of GSEA (gene set enrichment analysis) for RNA-Seq data
#	Yunshun Chen, Gordon Smyth
#	Created 20 Oct 2014.  Last modified 29 Dec 2017.
{
	if(is.null(design)) {
		design <- y$design
		if(is.null(design)) design <- model.matrix(~y$samples$group)
	}
	y <- .zscoreDGE(y=y, design=design, contrast=contrast)
	romer(y=y, index=index, design=design, contrast, ...)
}

.zscoreDGE <- function(y, design=NULL, contrast=ncol(design))
#	Calculate NB z-scores given a DGEList object.
#	Yunshun Chen, Gordon Smyth
#	Created 27 May 2015.
#	Last modified 1 May 2018.
{
#	Check for all zero counts
	allzero <- rowSums(y$counts>1e-8)==0
	if(any(allzero)) warning(sum(allzero),"rows with all zero counts")

#	Check dispersion estimates in y
	dispersion <- getDispersion(y)
	if(is.null(dispersion)) stop("Dispersion estimate not found. Please estimate the dispersion(s) before you proceed.")

#	Make default design matrix from group factor
	if(is.null(design)) design <- y$design
	if(is.null(design)) {
		if(nlevels(y$samples$group)<2) stop("design not supplied and samples all belong to the same group")
		design <- model.matrix(~y$samples$group)
		rownames(design) <- colnames(y)
	}
	nbeta <- ncol(design)
	if(nbeta < 2) stop("design matrix must have at least two columns")

#	contrast could be a coef name
	if(is.character(contrast)) {
		if(length(contrast)>1) stop("contrast should specify only one column of design")
		contrast <- which(contrast==colnames(design))
		if(!length(contrast)) stop("contrast doesn't match any column of design")
	}

#	Construct null hypothesis design matrix
	if(length(contrast) == 1) {
		design0 <- design[,-contrast,drop=FALSE]
	} else {
		design <- contrastAsCoef(design,contrast=contrast,first=FALSE)$design
		design0 <- design[,-nbeta,drop=FALSE]
	}

#	Null hypothesis fit
	fit.null <- glmFit(y, design0, prior.count=0)

#	Quantile residuals from null fit
#	Applying a floor to mu avoids infinite z-scores when mu=0
	y <- zscoreNBinom(y$counts, mu=pmax(fit.null$fitted.values,1e-17), size=1/dispersion)
	y
}
