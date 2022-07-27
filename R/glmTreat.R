glmTreat <- function(glmfit, coef=ncol(glmfit$design), contrast=NULL, lfc=log2(1.2), null="interval")
#	Likelihood ratio test or quasi-likelihood F-test with a threshold
#	Yunshun Chen and Gordon Smyth
#	Created on 05 May 2014. Last modified on 25 March 2018.
{
	if(lfc < 0) stop("lfc has to be non-negative")

#	Check if glmfit is from glmFit() or glmQLFit()
	isLRT <- is.null(glmfit$df.prior)
	Fit <- ifelse(isLRT, glmFit, glmQLFit)

#	Switch to glmLRT() or glmQLFTest() if lfc is zero
	if(lfc==0) {
		fun <- ifelse(isLRT, "glmLRT", "glmQLFTest")
		cat( paste0("Zero log2-FC threshold detected. Switch to ", fun, "() instead."), "\n" )
		return( do.call(fun, args=list(glmfit, coef, contrast)) )
	}

#	Check if there is any log-FC shrinkage
	shrunk <- glmfit$prior.count!=0
	
#	Check glmfit
	if(!is(glmfit,"DGEGLM")) {
		if(is(glmfit,"DGEList") && is(coef,"DGEGLM")) {
			stop("First argument is no longer required. Rerun with just the glmfit and coef/contrast arguments.")
		}
		stop("glmfit must be an DGEGLM object (usually produced by glmFit or glmQLFit).")
	}
	if(is.null(glmfit$AveLogCPM)) glmfit$AveLogCPM <- aveLogCPM(glmfit)
	ngenes <- nrow(glmfit)

#	Check design matrix
	design <- as.matrix(glmfit$design)
	nbeta <- ncol(design)
	if(nbeta < 2) stop("Need at least two columns for design, usually the first is the intercept column")
	coef.names <- colnames(design)

#	Evaluate logFC for coef to be tested
#	Note that contrast takes precedence over coef: if contrast is given
#	then reform design matrix so that contrast of interest is the first column
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
		unshrunk.logFC <- glmfit$coefficients[, coef, drop=FALSE]/log(2)
		if(shrunk) {
			logFC <- as.vector(unshrunk.logFC)
			unshrunk.logFC <- glmfit$unshrunk.coefficients[, coef, drop=FALSE]/log(2)
		}
	} else {
#		contrast should be a vector or a matrix with one column
		contrast <- as.matrix(contrast)
		if(ncol(contrast) > 1L) {
			warning("Found more than one contrast vector in `contrast`. Using first column only.")
			contrast <- contrast[,1,drop=FALSE]
		}
		reform <- contrastAsCoef(design, contrast=contrast, first=TRUE)
		coef <- 1
		unshrunk.logFC <- drop((glmfit$coefficients %*% contrast)/log(2))
		if(shrunk) {
			logFC <- as.vector(unshrunk.logFC)
			unshrunk.logFC <- drop((glmfit$unshrunk.coefficients %*% contrast)/log(2))
		}
		contrast <- drop(contrast)
		i <- contrast!=0
		coef.name <- paste(paste(contrast[i],coef.names[i],sep="*"),collapse=" ")
		design <- reform$design
	}
	unshrunk.logFC <- as.vector(unshrunk.logFC)

#	Null design matrix
	design0 <- design[, -coef, drop=FALSE]

#	Offset adjustment
	offset.old <- makeCompressedMatrix(glmfit$offset, dim(glmfit$counts), byrow=TRUE)
	offset.adj <- makeCompressedMatrix(lfc*log(2) * design[, coef], dim(glmfit$counts), byrow=TRUE)

#	Test statistics at beta_0 = tau
	offset.new <- offset.old + offset.adj
	fit0 <- Fit(glmfit$counts, design=design0, offset=offset.new, weights=glmfit$weights, dispersion=glmfit$dispersion, prior.count=0)
	fit1 <- Fit(glmfit$counts, design=design, offset=offset.new, weights=glmfit$weights, dispersion=glmfit$dispersion, prior.count=0)
	z.left <- sqrt( pmax(0, fit0$deviance - fit1$deviance) )

#	Test statistics at beta_0 = -tau
	offset.new <- offset.old - offset.adj
	fit0 <- Fit(glmfit$counts, design=design0, offset=offset.new, weights=glmfit$weights, dispersion=glmfit$dispersion, prior.count=0)
	fit1 <- Fit(glmfit$counts, design=design, offset=offset.new, weights=glmfit$weights, dispersion=glmfit$dispersion, prior.count=0)
	z.right <- sqrt( pmax(0, fit0$deviance - fit1$deviance) )

#	Make sure z.left < z.right
	i <- z.left > z.right
	if(any(i)) {
		tmp <- z.left[i]
		z.left[i] <- z.right[i]
		z.right[i] <- tmp
	}

#	Convert t to z under the QL pipeline
	if(!isLRT){
		df.total <- glmfit$df.prior + glmfit$df.residual.zeros
		max.df.residual <- ncol(glmfit$counts)-ncol(glmfit$design)
		df.total <- pmin(df.total, nrow(glmfit)*max.df.residual)
		z.left <- limma::zscoreT(z.left/sqrt(glmfit$var.post), df=df.total)
		z.right <- limma::zscoreT(z.right/sqrt(glmfit$var.post), df=df.total)
	}

	within <- abs(unshrunk.logFC) <= lfc
	sgn <- 2*within - 1
	z.left <- z.left*sgn

	null <- match.arg(null, c("interval", "worst.case"))
	if(null=="interval") {
#		Interval threshold
		c <- 1.470402
		j <- z.right + z.left > c
		p.value <- rep_len(1L, ngenes)
		j <- j[!is.na(j)]
		p.value[j] <- .integratepnorm(-z.right[j], -z.right[j] + c) + .integratepnorm(z.left[j] - c, z.left[j])
		p.value[!j] <- 2*.integratepnorm(-z.right[!j], z.left[!j])
	} else {
		p.value <- pnorm(-z.right) + pnorm(z.left)
	}
	
#	Ensure it is not more significant than chisquare test with Poisson variance		
	if(!isLRT){
		i <- rowSums(glmfit$var.post * (1 + glmfit$fitted.values * glmfit$dispersion) < 1) > 0L
		if(any(i)) {
			pois.fit <- glmfit[i,]
			pois.fit <- glmFit(pois.fit$counts, design=pois.fit$design, offset=pois.fit$offset, weights=pois.fit$weights, start=pois.fit$unshrunk.coefficients, dispersion=0)
			pois.res <- Recall(pois.fit, coef=coef, contrast=contrast, lfc=lfc)
			p.value[i] <- pmax(p.value[i], pois.res$table$PValue)
		}
	}

#	Table output
	tab <- data.frame(
		logCPM=glmfit$AveLogCPM,
		PValue=p.value,
		row.names=rownames(glmfit)
	)
	if(shrunk) {
		tab <- cbind(logFC=logFC, unshrunk.logFC=unshrunk.logFC, tab)
	} else {
		tab <- cbind(logFC=unshrunk.logFC, tab)
	}

	glmfit$lfc <- lfc
	glmfit$counts <- NULL
	glmfit$table <- tab
	glmfit$comparison <- coef.name
	new("DGELRT",unclass(glmfit))
}


.integratepnorm <- function(a, b) ifelse(a==b, pnorm(a), ( b*pnorm(b)+dnorm(b) - (a*pnorm(a)+dnorm(a)) )/(b-a) )
