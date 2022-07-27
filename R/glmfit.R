#  FIT GENERALIZED LINEAR MODELS

glmFit <- function(y, ...)
UseMethod("glmFit")

glmFit.DGEList <- function(y, design=NULL, dispersion=NULL, prior.count=0.125, start=NULL, ...)
#	Created 11 May 2011.  Last modified 21 Nov 2018.
{
#	The design matrix defaults to the oneway layout defined by y$samples$group.
#	If there is only one group, then the design matrix is left NULL so that a
#	matrix with a single intercept column will be set later by glmFit.default.
	if(is.null(design)) {
		design <- y$design
		if(is.null(design)) {
			group <- droplevels(as.factor(y$samples$group))
			if(nlevels(group) > 1L) design <- model.matrix(~y$samples$group)
		}
	}
	if(is.null(dispersion)) dispersion <- getDispersion(y)
	if(is.null(dispersion)) stop("No dispersion values found in DGEList object.")
	offset <- getOffset(y)
	if(is.null(y$AveLogCPM)) y$AveLogCPM <- aveLogCPM(y)

	fit <- glmFit(y=y$counts,design=design,dispersion=dispersion,offset=offset,lib.size=NULL,weights=y$weights,prior.count=prior.count,start=start,...)
	fit$samples <- y$samples
	fit$genes <- y$genes
	fit$prior.df <- y$prior.df
	fit$AveLogCPM <- y$AveLogCPM
	new("DGEGLM",fit)
}

glmFit.SummarizedExperiment <- function(y, design=NULL, dispersion=NULL, prior.count=0.125, start=NULL, ...)
#	Created 19 March 2020.  Last modified 19 March 2020.
{
	y <- SE2DGEList(y)
	glmFit.DGEList(y, design=design, dispersion=dispersion, prior.count=prior.count, start=start, ...)
}

glmFit.default <- function(y, design=NULL, dispersion=NULL, offset=NULL, lib.size=NULL, weights=NULL, prior.count=0.125, start=NULL, ...)
#	Fit negative binomial generalized linear model for each transcript
#	to a series of digital expression libraries
#	Davis McCarthy, Gordon Smyth, Yunshun Chen, Aaron Lun
#	Created 17 August 2010. Last modified 7 Aug 2019.
{
#	Check y
	y <- as.matrix(y)
	if(mode(y) != "numeric") stop("y is not a numeric matrix")
	ntag <- nrow(y)
	nlib <- ncol(y)

#	Check design
	if(is.null(design)) {
		design <- matrix(1,nlib,1)
		rownames(design) <- colnames(y)
		colnames(design) <- "Intercept"
	} else {
		design <- as.matrix(design)
		if(nrow(design) != nlib) stop("nrow(design) disagrees with ncol(y)")
		ne <- nonEstimable(design)
		if(!is.null(ne)) stop(paste("Design matrix not of full rank.  The following coefficients not estimable:\n", paste(ne, collapse = " ")))
	}

#	Check dispersion
	if(is.null(dispersion)) stop("No dispersion values provided.")
	if(anyNA(dispersion)) stop("NA dispersions not allowed")
	if(!is.numeric(dispersion)) stop("dispersion must be numeric")
	if(is.null(dim(dispersion))) {
		if( !any(length(dispersion)==c(1L,ntag)) ) stop("dispersion has wrong length. As a vector, it should agree with nrow(y)")
	} else {
		if( !all(dim(dispersion)==dim(y)) ) stop("Dimensions of dispersion don't agree with dimensions of y")
	}
	dispersion.mat <- .compressDispersions(y, dispersion)

#	Check offset
	if(!is.null(offset)) {
		if(!is.numeric(offset)) stop("offset must be numeric")
		if(is.null(dim(offset))) {
			if( !any(length(offset)==c(1L,nlib)) ) stop("offset has wrong length. As a vector, it should agree with ncol(y)")
		} else {
			if( !all(dim(offset)==dim(y)) ) stop("Dimensions of offset don't agree with dimensions of y")
		}
	}

#	Check lib.size
	if(!is.null(lib.size)) {
		if(!is.numeric(lib.size)) stop("lib.size must be numeric")
		if( !any(length(lib.size)==c(1L,nlib)) ) stop("lib.size has wrong length, should agree with ncol(y)")
	}

#	Comsolidate lib.size and offset into a compressed matrix
	offset <- .compressOffsets(y=y, lib.size=lib.size, offset=offset)

#	weights are checked in lower-level functions

#	Fit the tagwise glms
#	If the design is equivalent to a oneway layout, use a shortcut algorithm
	group <- designAsFactor(design)
	if(nlevels(group)==ncol(design)) {
		fit <- mglmOneWay(y,design=design,group=group,dispersion=dispersion.mat,offset=offset,weights=weights,coef.start=start)
		fit$deviance <- nbinomDeviance(y=y,mean=fit$fitted.values,dispersion=dispersion.mat,weights=weights)
		fit$method <- "oneway"
	} else {
		fit <- mglmLevenberg(y,design=design,dispersion=dispersion.mat,offset=offset,weights=weights,coef.start=start,maxit=250)
		fit$method <- "levenberg"
	}

#	Prepare output
	fit$counts <- y
	if(prior.count>0) {
		fit$unshrunk.coefficients <- fit$coefficients
		colnames(fit$unshrunk.coefficients) <- colnames(design)
		rownames(fit$unshrunk.coefficients) <- rownames(y)
		fit$coefficients <- predFC(y,design,offset=offset,dispersion=dispersion.mat,prior.count=prior.count,weights=weights,...)*log(2)
	}
	colnames(fit$coefficients) <- colnames(design)
	rownames(fit$coefficients) <- rownames(y)
	dimnames(fit$fitted.values) <- dimnames(y)
#	FIXME: we are not allowing missing values, so df.residual must be same for all tags
	fit$df.residual <- rep(nlib-ncol(design),ntag)
	fit$design <- design
	fit$offset <- offset
	fit$dispersion <- dispersion
	fit$weights <- weights
	fit$prior.count <- prior.count
	new("DGEGLM",fit)
}


glmLRT <- function(glmfit,coef=ncol(glmfit$design),contrast=NULL)
#	Tagwise likelihood ratio tests for DGEGLM
#	Gordon Smyth, Davis McCarthy and Yunshun Chen.
#	Created 1 July 2010.  Last modified 9 June 2020.
{
#	Check glmfit
	if(!is(glmfit,"DGEGLM")) {
		if(is(glmfit,"DGEList") && is(coef,"DGEGLM")) {
			stop("First argument is no longer required. Rerun with just the glmfit and coef/contrast arguments.")
		}
		stop("glmfit must be an DGEGLM object (usually produced by glmFit).")
	}
	if(is.null(glmfit$AveLogCPM)) glmfit$AveLogCPM <- aveLogCPM(glmfit)
	nlibs <- ncol(glmfit)
	
#	Check design matrix
	design <- as.matrix(glmfit$design)
	nbeta <- ncol(design)
	if(nbeta < 2) stop("Need at least two columns for design, usually the first is the intercept column")
	coef.names <- colnames(design)

#	Evaluate logFC for coef to be tested
#	Note that contrast takes precedence over coef: if contrast is given
#	then reform design matrix so that contrast of interest is last column.
	if(is.null(contrast)) {
		if(length(coef) > 1) coef <- unique(coef)
		if(is.character(coef)) {
			check.coef <- coef %in% colnames(design)
			if(any(!check.coef)) stop("One or more named coef arguments do not match a column of the design matrix.")
			coef.name <- coef
			coef <- match(coef, colnames(design))
		}
		else
			coef.name <- coef.names[coef]
		logFC <- glmfit$coefficients[,coef,drop=FALSE]/log(2)
	} else {
		contrast <- as.matrix(contrast)
		if(nrow(contrast) != ncol(glmfit$coefficients)) stop("contrast vector of wrong length, should be equal to number of coefficients in the linear model.")
		qrc <- qr(contrast)
		ncontrasts <- qrc$rank
		if(ncontrasts==0) stop("contrasts are all zero")
		coef <- 1:ncontrasts
		logFC <- (glmfit$coefficients %*% contrast)/log(2)
		if(ncontrasts>1) {
			coef.name <- paste("LR test on",ncontrasts,"degrees of freedom")
		} else {
			contrast <- drop(contrast)
			i <- contrast!=0
			coef.name <- paste(paste(contrast[i],coef.names[i],sep="*"),collapse=" ")
		}
		Dvec <- rep_len(1,nlibs)
		Dvec[coef] <- diag(qrc$qr)[coef]
		Q <- qr.Q(qrc,complete=TRUE,Dvec=Dvec)
		design <- design %*% Q
	}
	if(length(coef)==1) logFC <- drop(logFC)

#	Null design matrix
	design0 <- design[,-coef,drop=FALSE]

#	Null fit
	fit.null <- glmFit(glmfit$counts,design=design0,offset=glmfit$offset,weights=glmfit$weights,dispersion=glmfit$dispersion,prior.count=0)

#	Likelihood ratio statistic
	LR <- fit.null$deviance - glmfit$deviance
	df.test <- fit.null$df.residual - glmfit$df.residual
	LRT.pvalue <-  pchisq(LR, df=df.test, lower.tail = FALSE, log.p = FALSE)

	rn <- rownames(glmfit)
	if(is.null(rn))
		rn <- 1:nrow(glmfit)
	else
		rn <- make.unique(rn)
	tab <- data.frame(
		logFC=logFC,
		logCPM=glmfit$AveLogCPM,
		LR=LR,
		PValue=LRT.pvalue,
		row.names=rn
	)
	glmfit$counts <- NULL
	glmfit$table <- tab 
	glmfit$comparison <- coef.name
	glmfit$df.test <- df.test
	new("DGELRT",unclass(glmfit))
}

