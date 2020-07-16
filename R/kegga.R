kegga.DGELRT <- function(de, geneid = rownames(de), FDR = 0.05, trend = FALSE, ...)
#	KEGG analysis of DE genes from linear model fit
#	Gordon Smyth
#	Created 4 June 2015.  Last modified 27 May 2019.
{
#	Avoid argument collision with default method
	dots <- names(list(...))
	if("universe" %in% dots) stop("kegga.DGELRT defines its own universe",call.=FALSE)
	if((!is.logical(trend) || trend) && "covariate" %in% dots) stop("kegga.DGELRT defines it own covariate",call.=FALSE)
	ngenes <- nrow(de)

#	Check geneid
#	Can be either a vector of gene IDs or an annotation column name
	geneid <- as.character(geneid)
	if(length(geneid) == ngenes) {
		universe <- geneid
	} else {
		if(length(geneid) == 1L) {
			universe <- de$genes[[geneid]]
			if(is.null(universe)) stop("Column ",geneid," not found in de$genes")
		} else
			stop("geneid of incorrect length")
	}

#	Check trend
#	Can be logical, or a numeric vector of covariate values, or the name of the column containing the covariate values
	if(is.logical(trend)) {
		if(trend) {
			covariate <- de$table$logCPM
			if(is.null(covariate)) stop("logCPM not found in fit object")
		}
	} else {
		if(is.numeric(trend)) {
			if(length(trend) != ngenes) stop("If trend is numeric, then length must equal nrow(de)")
			covariate <- trend
			trend <- TRUE
		} else {
			if(is.character(trend)) {
				if(length(trend) != 1L) stop("If trend is character, then length must be 1")
				covariate <- de$genes[[trend]]
				if(is.null(covariate)) stop("Column ",trend," not found in de$genes")
				trend <- TRUE
			} else
				stop("trend is neither logical, numeric nor character")
		}
	}

#	Check FDR
	if(!is.numeric(FDR) | length(FDR) != 1) stop("FDR must be numeric and of length 1.")
	if(FDR < 0 | FDR > 1) stop("FDR should be between 0 and 1.")

#	If no DE genes, return data.frame with 0 rows
	sig <- (p.adjust(de$table$PValue, method = "BH") < FDR)
	if(sum(sig)==0L) {
		message("No DE genes")
		return(data.frame())
	}

#	Get up and down DE genes
	if(is.null(de$df.test)) de$df.test <- 1
	if(de$df.test[1] > 1) {
		DEGenes <- universe[sig]
	} else {
		EG.DE.UP <- universe[sig & de$table$logFC > 0]
		EG.DE.DN <- universe[sig & de$table$logFC < 0]
		DEGenes <- list(Up=EG.DE.UP, Down=EG.DE.DN)
	}

	if(trend)
		kegga(de=DEGenes, universe = universe, covariate=covariate, ...)
	else
		kegga(de=DEGenes, universe = universe, ...)
}


kegga.DGEExact <- kegga.DGELRT
