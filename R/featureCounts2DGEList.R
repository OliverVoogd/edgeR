featureCounts2DGEList <- function(x)
#	Convert featureCounts output to DGEList
#	Gordon Smyth
#	Created 7 Jan 2021. Last modified 26 Jan 2021.
{
#	Check for featureCounts output
	if(!is.list(x) || !all(c("counts","annotation","targets","stat") %in% names(x)))
		stop("x should be a list containing components `counts`, `annotation`, `targets` and `stat`")

#	Get and simplify sample names
	SampleNames <- colnames(x$counts)
	if(is.null(SampleNames)) SampleNames <- x$targets
	if(is.null(SampleNames)) SampleNames <- seq_len(ncol(x$counts))
	SampleNames <- removeExt(removeExt(SampleNames))
	colnames(x$counts) <- SampleNames

#	Mapping statistics	
	Stat <- t(as.matrix(x$stat[,-1,drop=FALSE]))
	rownames(Stat) <- SampleNames
	colnames(Stat) <- x$stat[,1]
	PropAssigned <- Stat[,1] / rowSums(Stat)
	Samples <- data.frame(PropAssigned=PropAssigned)

#	GeneID must be character
	x$annotation$GeneID <- as.character(x$annotation$GeneID)
	if(!identical(row.names(x$counts),x$annotation$GeneID)) warning("row.names of counts are not the same as GeneIDs")

#	Detect use of meta-features
	DupRowNames <- as.logical(anyDuplicated(rownames(x$counts)))
	MultipleFeatures <- as.logical(length(grep(";",x$annotation$Start)))
	Junctions <- !is.null(x$counts_junction)
	if(MultipleFeatures && DupRowNames) stop("Meta-feature counts should have unique row.names")
	if(MultipleFeatures && Junctions) stop("Junction counts incompatible with meta-feature counts")

##	Meta-features (genes)

	if(MultipleFeatures) {
#		Collapse Chr, Start, End and Strand for meta-feature
		firstfun <- function(z) z[1]
		lastfun <- function(z) z[length(z)]
		z <- strsplit(x$annotation$Chr,split=";")
		z1 <- x$annotation$Chr <- vapply(z,firstfun,"")
		z2 <- x$annotation$Chr <- vapply(z,lastfun,"")
		i <- z1 != z2
		z1[i] <- paste(z1[i],z2[i],sep=";")
		x$annotation$Chr <- z1
		z <- strsplit(as.character(x$annotation$Start),split=";")
		x$annotation$Start <- as.integer(vapply(z,firstfun,""))
		z <- strsplit(as.character(x$annotation$End),split=";")
		x$annotation$End <- as.integer(vapply(z,lastfun,""))
		z <- strsplit(x$annotation$Strand,split=";")
		x$annotation$Strand <- vapply(z,firstfun,"")

		x$annotation$GeneID <- NULL
		return(DGEList(counts=x$counts,genes=x$annotation,samples=Samples))
	}

##	Features

#	Ensure classes for features
	x$annotation$Start <- as.integer(x$annotation$Start)
	x$annotation$End <- as.integer(x$annotation$End)
	x$annotation$Strand <- as.character(x$annotation$Strand)
	x$annotation$Length <- as.integer(x$annotation$Length)

##	Simple features

	if(!DupRowNames && !Junctions) {
		x$annotation$GeneID <- NULL
		return(DGEList(counts=x$counts,genes=x$annotation,samples=Samples))
	}

##	Features within meta-features (exons)

#	Expand the row names to include exon number (counting in genomic order)
	o <- order(x$annotation$GeneID,x$annotation$Chr,x$annotation$Start)
	Counts <- x$counts[o,]
	Genes <- x$annotation[o,]
	FirstExon <- which(!duplicated(Genes$GeneID))
	LastExon <- c(FirstExon[-1]-1L,nrow(Counts))
	RowNum <- seq_len(nrow(Counts))
	NExons <- RowNum[LastExon]-RowNum[FirstExon]+1L
	ExonNum <- RowNum-RowNum[rep.int(FirstExon,times=NExons)]+1L
	Genes$Exon <- ExonNum
	rownames(Counts) <- paste0(Genes$GeneID,".e",ExonNum)
	row.names(Genes) <- rownames(Counts)

#	Sort in genomic order	
	o <- order(Genes$Chr,Genes$Start)
	Counts <- Counts[o,]
	Genes <- Genes[o,]

	if(!Junctions) return (DGEList(counts=Counts,genes=Genes,samples=Samples))

##	Exons and junctions

#	JCounts <- as.matrix(x$counts_junction[,-seq_len(8)])
#	colnames(JCounts) <- SampleNames
#	JGenes <- x$counts_junction[,c("PrimaryGene","Site1_chr","Site1_location","Site2_location","Site1_strand")]
#	names(JGenes) <- c("GeneID","Chr","Start","End","Strand")
#	JGenes$GeneID <- as.character(JGenes$GeneID)
#	JGenes$Length <- JGenes$End-JGenes$Start+1L
#	JGenes$Exon <- rep_len(0L,nrow(JGenes))

	JCounts <- x$counts_junction
    JCounts <- JCounts[!is.na(JCounts$PrimaryGene),]
    JGenes <- JCounts[,c("PrimaryGene","Site1_chr","Site1_location","Site2_location","Site1_strand")]
    names(JGenes) <- c("GeneID","Chr","Start","End","Strand")
	JGenes$GeneID <- as.character(JGenes$GeneID)
	JGenes$Start <- as.integer(JGenes$Start)
	JGenes$End <- as.integer(JGenes$End)
	JGenes$Strand <- as.character(JGenes$Strand)
    JGenes$Length <- 1L
    JCounts <- as.matrix(JCounts[,-seq_len(8)])
    o <- order(JGenes$GeneID,JGenes$Chr,JGenes$Start)
    JCounts <- JCounts[o,]
    JGenes <- JGenes[o,]
    RowNames <- JGenes$GeneID
    N <- length(RowNames)
    First <- which(!duplicated(RowNames))
    Last <- c(First[-1]-1L,N)
    NExons <- Last-First+1L
    One <- rep_len(1L,N)
    One[First[-1]] <- 1L-NExons[-length(NExons)]
    ExonNumber <- cumsum(One)
    RowNames <- paste0(RowNames,".j",ExonNumber)
    JGenes$Exon <- ExonNumber
    row.names(JGenes) <- RowNames
    rownames(JCounts) <- RowNames

	DGEList(counts=rbind(Counts,JCounts),genes=rbind(Genes,JGenes),samples=Samples)
}
