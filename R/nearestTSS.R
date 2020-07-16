nearestTSS <- function(chr,locus,species="Hs")
#	Find nearest gene transcriptional start sites from the appropriate organism db package
#	Gordon Smyth
#	Created 3 Jan 2018.  Last modified 21 Dec 2019.
{
#	Check input
	chr <- as.character(chr)
	locus <- as.integer(locus)
	n <- length(chr)
	if(length(locus) == 1L) 
		locus <- rep_len(locus,n)
	else
		if(length(locus) != n) stop("Length of locus doesn't agree with length of chr")

#	Replace NAs with empty string
	if(anyNA(chr)) {
		chr[is.na(chr)] <- ""
	}
	if(anyNA(locus)) {
		chr[is.na(locus)] <- ""
		locus[is.na(locus)] <- 0L
	}

#	Get access to required annotation functions
	suppressPackageStartupMessages(OK <- requireNamespace("AnnotationDbi",quietly=TRUE))
	if(!OK) stop("AnnotationDbi package required but is not installed (or can't be loaded)")

#	Load appropriate organism package
	orgPkg <- paste0("org.",species,".eg.db")
	suppressPackageStartupMessages(OK <- requireNamespace(orgPkg,quietly=TRUE))
	if(!OK) stop(orgPkg," package required but is not installed (or can't be loaded)")

#	Get gene start positions
	obj <- paste0("org.",species,".egCHRLOC")
	egCHRLOC <- tryCatch(getFromNamespace(obj,orgPkg), error=function(e) FALSE)
	if(is.logical(egCHRLOC)) stop("Can't find egCHRLOC gene location mappings in package ",orgPkg)
	EGLOC <- AnnotationDbi::toTable(egCHRLOC)

#	Get first gene end position
	obj <- paste0("org.",species,".egCHRLOCEND")
	egCHRLOCEND <- tryCatch(getFromNamespace(obj,orgPkg), error=function(e) FALSE)
	if(is.logical(egCHRLOCEND)) stop("Can't find egCHRLOCEND gene end mappings in package ",orgPkg)
	EGEND <- AnnotationDbi::toTable(egCHRLOCEND)
	EGLOC$end_location <- EGEND$end_location

#	Get Symbols
	obj <- paste0("org.",species,".egSYMBOL")
	egSYMBOL <- tryCatch(getFromNamespace(obj,orgPkg), error=function(e) FALSE)
	if(is.logical(egSYMBOL)) stop("Can't find egSYMBOL gene symbol mappings in package ",orgPkg)
	EGSym <- AnnotationDbi::toTable(egSYMBOL)
	m <- match(EGLOC$gene_id,EGSym$gene_id)
	EGLOC$symbol <- EGSym[m,2]

#	Get strand, width and TSS
	EGLOC$neg <- (EGLOC$start_location < 0L)
	EGLOC$width <- EGLOC$end_location - EGLOC$start_location
	EGLOC$tss <- EGLOC$start_location+1L
	EGLOC$tss[EGLOC$neg] <- EGLOC$end_location[EGLOC$neg]

#	Keep only positive loci
	EGLOC$width <- abs(EGLOC$width)
	EGLOC$tss <- abs(EGLOC$tss)
	EGLOC$start_location <- EGLOC$end_location <- NULL
	EGLOC$strand <- rep_len("+",nrow(EGLOC))
	EGLOC$strand[EGLOC$neg] <- "-"

#	Sort by genomic position
	o <- order(EGLOC$Chromosome,EGLOC$tss)
	EGLOC <- EGLOC[o,]
#	NREF <- nrow(EGLOC)
#	ChrFirst <- which(EGLOC$Chromosome[-NREF] != EGLOC$Chromosome[-1L])
#	ChrN <- c(1L,ChrFirst+1L)

#	Do chr values start with "chr"?
	if(length(grep("^chr",chr[1]))) EGLOC$Chromosome <- paste0("chr",EGLOC$Chromosome)

#	Prepare output
	n <- length(chr)
	ILocus <- rep_len(0L,n)

#	Cycle over chromosomes
	ChrNames <- unique(EGLOC$Chromosome)
	for (ChrA in ChrNames) {
		iref <- which(EGLOC$Chromosome==ChrA)
		iinc <- which(chr==ChrA)
		Which <- nearestReftoX(locus[iinc], EGLOC$tss[iref])
		ILocus[iinc] <- iref[Which]
	}

#	Reorder EGLOC to match input loci
	EGLOC$Chromosome <- NULL
	Out <- EGLOC[ILocus,,drop=FALSE]
	Out$distance <- Out$tss - locus[ILocus > 0L]
	Out$distance[Out$neg] <- -Out$distance[Out$neg]
	Out$neg <- NULL

#	If any incoming loci not found, return rows with NAs
	if(nrow(Out) < n) {
		ChrNA <- rep_len(NA_character_,n)
		IntNA <- rep_len(NA_integer_,n)
		Out2 <- data.frame(gene_id=ChrNA,symbol=ChrNA,width=IntNA,tss=IntNA,strand=ChrNA,distance=IntNA,stringsAsFactors=FALSE)
		Out2[ILocus > 0L,] <- Out
		return(Out2)
	} else {
		row.names(Out) <- 1:n
		return(Out)
	}
}
