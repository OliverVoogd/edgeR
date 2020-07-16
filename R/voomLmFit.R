voomLmFit <- function(
	counts, design=NULL, block=NULL, prior.weights=NULL,
	sample.weights=FALSE, var.design=NULL, var.group=NULL, 
	lib.size=NULL, normalize.method="none",
	span=0.5, plot=FALSE, save.plot=FALSE
)
#	limma+lmFit pipeline for counts taking into account of structural zeros
#	Creates an MArrayLM object for entry to eBayes() etc in the limma pipeline.
#	Depends on edgeR as well as limma
#	Gordon Smyth
#	Created 21 Jan 2020.  Last modified 10 Jun 2020.
{
	Block <- !is.null(block)
	PriorWeights <- !is.null(prior.weights)
	SampleWeights <- sample.weights || !is.null(var.design) || !is.null(var.group)

#	Can't specify prior weights and ask for sample weights to be estimated as well
	if(PriorWeights && SampleWeights) stop("Can't specify prior.weights and estimate sample weights")

#	Create output object
	out <- list()

#	Extract counts from known data objects
	if(is(counts,"SummarizedExperiment")) counts <- SE2DGEList(counts)
	if(is(counts,"DGEList")) {
		out$genes <- counts$genes
		out$targets <- counts$samples
		if(is.null(design) && diff(range(as.numeric(counts$sample$group)))>0) design <- model.matrix(~group,data=counts$samples)
		if(is.null(lib.size)) lib.size <- effectiveLibSizes(counts)
		counts <- counts$counts
	} else {
		if(is(counts,"eSet")) {
			if(!requireNamespace("Biobase",quietly=TRUE))
				stop("Biobase package required but is not installed (or can't be loaded)")
			if(length(Biobase::fData(counts))) out$genes <- Biobase::fData(counts)
			if(length(Biobase::pData(counts))) out$targets <- Biobase::pData(counts)
			counts <- get("counts",Biobase::assayData(counts))
		} else {
			counts <- as.matrix(counts)
		}
	}

#	Check counts
	n <- nrow(counts)
	if(n < 2L) stop("Need at least two genes to fit a mean-variance trend")
	m <- min(counts)
	if(is.na(m)) stop("NA counts not allowed")
	if(m < 0) stop("Negative counts not allowed")

#	Check design
	if(is.null(design)) {
		design <- matrix(1,ncol(counts),1)
		rownames(design) <- colnames(counts)
		colnames(design) <- "GrandMean"
	}

#	Check lib.size
	if(is.null(lib.size)) lib.size <- colSums(counts)

#	Expand prior.weights if necessary
	if(!is.null(prior.weights)) prior.weights <- asMatrixWeights(prior.weights,dim(counts))

#	log2-counts-per-million
	y <- t(log2(t(counts+0.5)/(lib.size+1)*1e6))

#	Microarray-style normalization
	y <- normalizeBetweenArrays(y,method=normalize.method)

#	Fit linear model
	fit <- lmFit(y,design,weights=prior.weights)

#	Find largest leverage value of design matrix
	if(is.null(fit$qr))
		h <- hat(design,intercept=FALSE)
	else
		h <- hat(fit$qr)
	MinGroupSize <- 1/max(h)

#	Identify fitted values that are exactly zero and should not contribute to the genewise variances
#	Note that a single zero is never a problem
	RowHasZero <- which(rowSums(counts==0) >= max(2,MinGroupSize))
	AnyZeroRows <- as.logical(length(RowHasZero))
	if(AnyZeroRows) {
		countsZero <- counts[RowHasZero,,drop=FALSE]
		PoissonFit <- glmFit(countsZero,design=design,lib.size=lib.size,dispersion=0,prior.count=0)
		IsZero <- (PoissonFit$fitted.values < 1e-4 & countsZero < 1e-4)
		RowHasExactZero <- which(rowSums(IsZero) > 0)
#		If any exact zero fits, then rerun the linear model for those rows with NAs
		if(length(RowHasExactZero)) {
			RowHasZero <- RowHasZero[RowHasExactZero]
			IsZero <- IsZero[RowHasExactZero,,drop=FALSE]
			yNAshort <- y[RowHasZero,,drop=FALSE]
			yNAshort[IsZero] <- NA
			fitNA <- suppressWarnings(lmFit(yNAshort,design,weights=prior.weights[RowHasZero,,drop=FALSE]))
			fit$df.residual[RowHasZero] <- fitNA$df.residual
			fit$sigma[RowHasZero] <- fitNA$sigma
#			If blocking or sample weights are present, then we will later on need a full length copy of y with NAs inserted
			if(Block || SampleWeights) {
				yNAfull <- y
				yNAfull[RowHasZero,] <- yNAshort
			}
		} else {
			AnyZeroRows <- FALSE
		}
	}

#	If no replication found, assume all weights are 1 and return fit already computed
	HasRep <- (fit$df.residual > 0L)
	NWithReps <- sum(HasRep)
	if(NWithReps < 2L) {
		if(NWithReps == 0L) warning("The experimental design has no replication. Setting weights to 1.")
		if(NWithReps == 1L) warning("Only one gene with any replication. Setting weights to 1.")
		fit$genes <- out$genes
		return(fit)
	}

#	Fit lowess trend to sqrt-standard-deviations by log-count-size
	Amean <- Amean2 <- rowMeans(y)
	if(AnyZeroRows) Amean2[RowHasZero] <- rowMeans(yNAshort,na.rm=TRUE)
	sx <- Amean2[HasRep]+mean(log2(lib.size+1))-log2(1e6)
	sy <- sqrt(fit$sigma[HasRep])
	if(AnyZeroRows)
		l <- weightedLowess(sx,sy,span=span,weights=fit$df.residual[HasRep],output.style="lowess")
	else
		l <- lowess(sx,sy,f=span)
	if(plot) {
		plot(sx,sy,xlab="log2( count size + 0.5 )",ylab="Sqrt( standard deviation )",pch=16,cex=0.25)
		title("voom: Mean-variance trend")
		lty <- ifelse(Block || SampleWeights,2,1)
		lines(l,col="red",lty=lty)
	}

#	Make interpolating rule
	f <- approxfun(l, rule=2, ties=list("ordered",mean))

#	Find individual quarter-root fitted counts
	if(fit$rank < ncol(design)) {
		j <- fit$pivot[1:fit$rank]
		fitted.values <- fit$coefficients[,j,drop=FALSE] %*% t(fit$design[,j,drop=FALSE])
	} else {
		fitted.values <- fit$coefficients %*% t(fit$design)
	}
	fitted.cpm <- 2^fitted.values
	fitted.count <- 1e-6 * t(t(fitted.cpm)*(lib.size+1))
	fitted.logcount <- log2(fitted.count)

#	Apply trend to individual observations to get voom weights
	w <- 1/f(fitted.logcount)^4
	dim(w) <- dim(fitted.logcount)

#	Add voom weights to prior weights
	if(PriorWeights)
		weights <- w * prior.weights
	else
		weights <- w

#	Estimate sample weights?
	if(SampleWeights) {
		if(AnyZeroRows) {
			sw <- arrayWeights(yNAfull,design,weights=weights,var.design=var.design,var.group=var.group)
		} else {
			sw <- arrayWeights(y,design,weights=weights,var.design=var.design,var.group=var.group)
		}
		message("First sample weights (min/max) ", paste(format(range(sw)),collapse="/") )
		if(Block) weights <- t(sw * t(weights))
	}

#	Estimate correlation?
	if(Block) {
		if(AnyZeroRows) {
			dc <- suppressWarnings(duplicateCorrelation(yNAfull,design,block=block,weights=weights))
		} else {
			dc <- suppressWarnings(duplicateCorrelation(y,design,block=block,weights=weights))
		}
		correlation <- dc$consensus.correlation
		if(is.na(correlation)) correlation <- 0
		message("First intra-block correlation  ",format(correlation))
	} else {
		correlation <- NULL
	}

#	Seond iteration to refine intra-block correlation or sample weights
	if(Block || SampleWeights) {
#		Rerun voom weights with new correlation and sample weights
		if(SampleWeights)
			weights <- asMatrixWeights(sw,dim(y))
		else
			weights <- prior.weights
		fit <- lmFit(y,design,block=block,correlation=correlation,weights=weights)
		if(AnyZeroRows) {
			fitNA <- suppressWarnings(lmFit(yNAshort,design,block=block,correlation=correlation,weights=weights[RowHasZero,,drop=FALSE]))
			fit$df.residual[RowHasZero] <- fitNA$df.residual
			fit$sigma[RowHasZero] <- fitNA$sigma
		}
		sy <- sqrt(fit$sigma[HasRep])
		if(AnyZeroRows)
			l <- weightedLowess(sx,sy,span=span,weights=fit$df.residual[HasRep],output.style="lowess")
		else
			l <- lowess(sx,sy,f=span)
		if(plot) {
			lines(l,col="red")
			legend("topright",lty=c(2,1),col="red",legend=c("First","Final"))
		}
		f <- approxfun(l, rule=2, ties=list("ordered",mean))
		if(fit$rank < ncol(design)) {
			j <- fit$pivot[1:fit$rank]
			fitted.values <- fit$coefficients[,j,drop=FALSE] %*% t(fit$design[,j,drop=FALSE])
		} else {
			fitted.values <- fit$coefficients %*% t(fit$design)
		}
		fitted.cpm <- 2^fitted.values
		fitted.count <- 1e-6 * t(t(fitted.cpm)*(lib.size+1))
		fitted.logcount <- log2(fitted.count)
		w <- 1/f(fitted.logcount)^4
		dim(w) <- dim(fitted.logcount)
		if(PriorWeights)
			weights <- w * prior.weights
		else
			weights <- w
		if(SampleWeights) {
			if(AnyZeroRows) {
				sw <- arrayWeights(yNAfull,design,weights=weights,var.design=var.design,var.group=var.group)
			} else {
				sw <- arrayWeights(y,design,weights=weights,var.design=var.design,var.group=var.group)
			}
			message("Final sample weights (min/max) ", paste(format(range(sw)),collapse="/") )
			weights <- t(sw * t(weights))
		}
		if(Block) {
			if(AnyZeroRows) {
				dc <- suppressWarnings(duplicateCorrelation(yNAfull,design,block=block,weights=weights))
			} else {
				dc <- suppressWarnings(duplicateCorrelation(y,design,block=block,weights=weights))
			}
			correlation <- dc$consensus.correlation
			if(is.na(correlation)) correlation <- 0
			message("Final intra-block correlation  ",format(correlation))
		}
	}

#	Final linear model fit with voom weights
	fit <- lmFit(y,design,block=block,correlation=correlation,weights=weights)
	if(is.null(fit$Amean)) fit$Amean <- Amean
	if(AnyZeroRows) {
		fitNA <- suppressWarnings(lmFit(yNAshort,design,block=block,correlation=correlation,weights=weights[RowHasZero,,drop=FALSE]))
		fit$df.residual[RowHasZero] <- fitNA$df.residual
		fit$sigma[RowHasZero] <- fitNA$sigma
	}

#	Output
	fit$genes <- out$genes
	fit$targets <- out$targets
	if(is.null(fit$targets)) {
		fit$targets <- data.frame(lib.size=lib.size)
		row.names(fit$targets) <- colnames(y)
	}
	if(SampleWeights) fit$targets$sample.weights <- sw
	if(save.plot) {
		fit$voom.xy <- list(x=sx,y=sy,xlab="log2( count size + 0.5 )",ylab="Sqrt( standard deviation )")
		fit$voom.line <- l
	}
	fit
}
