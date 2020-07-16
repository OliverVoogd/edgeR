spliceVariants <- function(y, geneID, dispersion=NULL, group=NULL, estimate.genewise.disp=TRUE, trace=FALSE)
# Identify genes with splice variants using a negative binomial model
# We assume that the data come in a matrix (possibly and/or a DGEList), counts summarized at exon level, with gene information available
# Davis McCarthy and Gordon Smyth
# Created 4 February 2011.  Last modified 3 Oct 2016. 
{
	if( is(y, "DGEList") ) {
		y.mat <- y$counts
		if( is.null(group) )
			group <- y$samples$group
		lib.size <- y$samples$lib.size * y$samples$norm.factors
	}
	else {
		y.mat <- as.matrix(y)
		if( is.null(group) )
			stop("y is a matrix and no group argument has been supplied. Please supply group argument.")
		lib.size <- colSums(y)
	}

	geneID <- as.vector(unlist(geneID))
#	Order genes by geneID: we need some way to reorganise the data---output cannot possibly be same dimension as input so this is a sensible way to organise things
	o <- order(geneID)
	geneID <- geneID[o]
	y.mat <- y.mat[o,]
	uniqIDs <- unique(geneID)
#	Organise the dispersion values to be used in the NB models
	if( is.null(dispersion) ) {
		if(estimate.genewise.disp) {
			dispersion <- estimateExonGenewiseDisp(y.mat, geneID, group)
			genewise.disp <- TRUE
			if(trace)
				cat("Computing genewise dispersions from exon counts.\n")
		}
		else {
			dispersion <- estimateCommonDisp(DGEList(counts=y.mat, group=group))
			genewise.disp <- FALSE
			if(trace)
				cat("Computing common dispersion across all exons.\n")
		}
	}
	else {
		if( length(dispersion)==length(uniqIDs) ) {
			if(is.null(names(dispersion)))
				stop("names(dispersion) is NULL. All names of dispersion must be unique geneID.\n")
			matches <- match(uniqIDs, names(dispersion))
			if( any(is.na(matches)) | any(duplicated(matches)))
			   stop("names(dispersion) of dispersion do not have a one-to-one mapping to unique geneID. All names of dispersion must be unique geneID.\n")
			dispersion <- dispersion[matches]
			genewise.disp <- TRUE
		}
		if( length(dispersion)==1 ) {
			genewise.disp <- FALSE
		}
	}
#	Can't get any results if there are no counts for an exon, so remove all-zero exons
	keep <- rowSums(y.mat) > 0
	exons <- y.mat[keep,]
	rownames(exons) <- geneID[keep]
	uniqIDs <- unique(geneID[keep])
	na.vec <- rep_len(NA_real_, length(uniqIDs))
	if(genewise.disp)
		dispersion <- dispersion[names(dispersion) %in% uniqIDs]
	else {
		dispersion <- rep_len(dispersion, length(uniqIDs))
		names(dispersion) <- uniqIDs
	}
#	We want to know how many exons each gene has
	dummy <- rowsum(rep_len(1, nrow(exons)), rownames(exons))
	nexons <- as.vector(dummy)
	names(nexons) <- rownames(dummy)
	mm <- match(uniqIDs,names(nexons))
	nexons <- nexons[mm]
	if(trace)
		cat("Max number exons: ",max(nexons),"\n")
#	Genes with the same number of exons have the same design matrix, allowing some parallelization of computations
	splicevars.out <- data.frame(logFC= na.vec, logCPM = na.vec, LR = na.vec, PValue = na.vec)
	rownames(splicevars.out) <- uniqIDs
	abundance <- na.vec
#	For loop iterates over number of exons for genes, starting at 2 (can't have splice variants if only one exon!)
	for(i.exons in sort(unique(nexons))) {
	#	Select the genes with the right number of exons for this iteration
		this.genes <- nexons==i.exons
		full.index <- rownames(exons) %in% uniqIDs[this.genes]
		if( any(this.genes) ) {
			gene.counts.mat <- matrix(t(exons[full.index,]), nrow=sum(this.genes), ncol=ncol(exons)*i.exons, byrow=TRUE)
			expanded.lib.size <- rep.int(lib.size, i.exons)
			if(i.exons==1) {
				abundance[this.genes] <- aveLogCPM(gene.counts.mat, lib.size=expanded.lib.size)
				splicevars.out$LR[this.genes] <- 0
				splicevars.out$PValue[this.genes] <- 1
			}
			else {
				exon.this <- factor(rep(1:i.exons, each=ncol(exons)))
				group.this <- as.factor(rep.int(group, i.exons))
			#	Define design matrices for this group of genes
				X.full <- model.matrix(~ exon.this + group.this + exon.this:group.this )
				X.null <- model.matrix(~ exon.this + group.this )
				coef <- (ncol(X.null)+1):ncol(X.full)
			#	Fit NB GLMs to these genes
				fit.this <- glmFit(gene.counts.mat, X.full, dispersion[this.genes], offset=0, prior.count=0)
				abundance[this.genes] <- aveLogCPM(gene.counts.mat, lib.size=expanded.lib.size)
				results.this <- glmLRT(fit.this, coef=coef)
				if(sum(this.genes)==1) {
					splicevars.out$LR[this.genes] <- results.this$table$LR[1]
					splicevars.out$PValue[this.genes] <- results.this$table$PValue[1]
				}
				else {
					splicevars.out$LR[this.genes] <- results.this$table$LR
					splicevars.out$PValue[this.genes] <- results.this$table$PValue
				}
			}
		}
	}
	splicevars.out$logCPM <- abundance
#	Create a list with the exons divided up neatly by geneID (a bit slow using in-built fn)
#	Not really necessary, so leave out for the time being
#	exon.list <- split(as.data.frame(exons), rownames(exons))
#	names(exon.list) <- uniqIDs
	if(!genewise.disp)
		dispersion <- dispersion[1]
	new("DGEExact",list(table=splicevars.out, comparison=NULL, genes=data.frame(GeneID=uniqIDs), dispersion=dispersion))
}



estimateExonGenewiseDisp <- function(y, geneID, group=NULL)
#	Function to estimate a common dispersion from exon count data
#	Created by Davis McCarthy, 29 July 2011.
#	Last modified 7 Aug 2019.
{
#	Check objects coming in
	if( is(y, "DGEList") ) {
		y.mat <- y$counts
		if( is.null(group) )
			group <- y$samples$group
	}
	else {
		y.mat <- as.matrix(y)
		if( is.null(group) )
			stop("y is a matrix and no group argument has been supplied. Please supply group argument\n")
	}
	geneID <- as.vector(unlist(geneID))
#	Cannot maintain order of the argument y, so order on geneID so that we have some sensible organisation
	o <- order(geneID)
	geneID <- geneID[o]
	y.mat <- y.mat[o,]
#	Sum counts from all exons for each gene to get gene-level counts and form DGEList object
	gene.counts <- rowsum(y.mat, geneID)
	genewise.disp <- rep_len(NA_real_, nrow(gene.counts))
	names(genewise.disp) <- rownames(gene.counts)
	gene.data <- DGEList(counts=gene.counts, group=group)
#	Cannot properly compute dispersion if there are no counts for the gene
	used <- rowSums(gene.data$counts) > 0
#	Need to first estimate a common dispersion
	gene.data <- estimateCommonDisp(gene.data[used,])
#	Next estimate tagwise dispersion for each gene, with trend. Default prop.used=2/3, grid.length=200
	gene.data <- estimateTagwiseDisp(gene.data, trend="movingave")
#	For those gene which have sufficient (>0) counts, assign estimated dispersion
	genewise.disp[used] <- gene.data$tagwise.dispersion
#	For those genes with zero counts, assign maximum estimated dispersion value
	genewise.disp[!used] <- max(gene.data$tagwise.dispersion)
	genewise.disp
}


plotExonUsage <- function(y, geneID, group=NULL, transform="none", counts.per.million=TRUE, legend.coords=NULL, ...)
#	Plots exon usage from a matrix, DGEList or list of exon counts
#	Created by Davis McCarthy, 30 July 2011.
#	Last modified 2 Aug 2011.
{
	if( is(y,"DGEList") ) {
		ind <- rownames(y$counts) %in% geneID
		counts <- y$counts
		if(counts.per.million)
			counts <- cpm(counts)
		exon.mat <- counts[ind,]
		group <- y$samples$group
	}
	else {
		if( is(y,"list") ) {
			exon.mat <- y[[geneID]]
			if(counts.per.million)
				stop("Counts per million cannot easily be computed when y is a list. Please use either a DGEList object or a matrix of counts.\n")
			if(is.null(group))
				stop("Group argument must be supplied.\n")
		}
		else {
			if(is.null(group))
				stop("Group argument must be supplied.\n")
			if(counts.per.million)
				y <- cpm(y)
			ind <- rownames(y) %in% geneID
			exon.mat <- y[ind,]
		}
	}
	transform <- match.arg(transform, c("none", "log2", "sqrt"))
	if(transform=="none") {
		if(counts.per.million)
			ylab <- "Counts per million"
		else
			ylab <- "Read counts"
	}
	if(transform=="log2") {
		exon.mat <- log2(exon.mat + 0.5)
		if(counts.per.million)
			ylab <- "log2( Counts per million )"
		else
			ylab <- "log2( Read counts )"
	}
	if(transform=="sqrt") {
		exon.mat <- sqrt(exon.mat)
		if(counts.per.million)
			ylab <- "sqrt( Counts per million )"
		else
			ylab <- "sqrt( Read counts )"
	}
	if(length(group)!=ncol(exon.mat))
		stop("Length of group vector does not match number of samples (columns of exon matrix).\n")
	cols <- 1:nlevels(group)
	plot(exon.mat[,1], type="b", col=cols[group[1]], panel.first=grid(), ylab=ylab, xlab="Exon", main=paste("GeneID: ",geneID), ylim=c(0, max(exon.mat)), ...)
	for( i in 2:ncol(exon.mat) )
		lines(exon.mat[,i], type="b", col=cols[group[i]], ...)
	if(is.null(legend.coords)) {
		legend.coords <- rep_len(NA_real_,2L)
		legend.coords[1] <- 0.8*nrow(exon.mat)
		legend.coords[2] <- 0.9*max(exon.mat)
	}
	legend(legend.coords[1], legend.coords[2], levels(group), col=cols, lwd=2)
	invisible(exon.mat)
}

