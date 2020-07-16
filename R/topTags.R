topTags <- function(object,n=10L,adjust.method="BH",sort.by="PValue",p.value=1) 
#	Summary table of the n most differentially expressed tags
#	Mark Robinson, Davis McCarthy, Gordon Smyth, Yunshun Chen.
#	Created September 2008.  Last modified 9 Oct 2017.
{
#	Check object
	if(is.null(object$table)) stop("Need to run exactTest or glmLRT first")
	if(is(object,"DGEExact")) test <- "exact" else test <- "glm"
	MultipleContrasts <- (test=="glm" && ncol(object$table) > 4)

#	Check n
	n <- min(n,nrow(object$table))
	if(n<1L) stop("No rows to output")

#	Check adjust.method
	FWER.methods <- c("holm", "hochberg", "hommel", "bonferroni")
	FDR.methods <- c("BH", "BY", "fdr")
	adjust.method <- match.arg(adjust.method,c(FWER.methods,FDR.methods,"none"))

#	Check sort.by
	sort.by <- match.arg(sort.by,c("none","p.value","PValue","logFC"))
	if(sort.by=="p.value") sort.by <- "PValue"

#	Absolute log fold change
	if(MultipleContrasts) {
		if(sort.by=="logFC") warning("Two or more logFC columns in DGELRT object. First logFC column used to rank by logFC.")
		alfc <- abs(object$table[,1])
	} else {
		alfc <- abs(object$table$logFC)
	}

#	Choose top genes
	o <- switch(sort.by,
		"logFC" = order(alfc,decreasing=TRUE),
		"PValue" = order(object$table$PValue,-alfc),
		"none" = 1:nrow(object$table)
	)
	tab <- object$table[o,]

#	Add adjusted p-values if appropriate
	adj.p.val <- p.adjust(object$table$PValue,method=adjust.method)
	if(adjust.method != "none") {
		if(adjust.method %in% FWER.methods) adjustment <- "FWER"
		if(adjust.method %in% FDR.methods) adjustment <- "FDR"
		tab[[adjustment]] <- adj.p.val[o]
	}

#	Add gene annotation if appropriate
	if(!is.null(object$genes)){
		if(is.null(dim(object$genes))) object$genes <- data.frame(ID=object$genes,stringsAsFactors=FALSE)
		rn <- row.names(tab)
		tab <- cbind(object$genes[o,,drop=FALSE], tab)
		row.names(tab) <- rn
	}
	
#	Thin out fit p.value threshold
	if(p.value < 1) {
		sig <- adj.p.val[o] <= p.value
		sig[is.na(sig)] <- FALSE
		tab <- tab[sig,]
	}

#	Enough rows left?
	if(nrow(tab) < n) n <- nrow(tab)
	if(n < 1L) return(data.frame())
		
#	Output object
	new("TopTags",list(
		table=tab[1:n,],
		adjust.method=adjust.method,
		comparison=as.character(object$comparison),
		test=test
	))
}

