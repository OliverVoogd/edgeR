readBismark2DGE <- function(files,sample.names=NULL,readr=TRUE,verbose=TRUE)
#	Read Bismark coverage files and create DGEList
#
#	It is assumed that genomic loci can be represented as integers, so
#	the largest locus position must be less than about 2*10^9.
#	The number of chromosomes times the largest locus position must be
#	less than 10^16.
#
#	Gordon Smyth
#	Created 30 May 2018. Last modified 26 Feb 2020.
{
	files <- as.character(files)
	nsamples <- length(files)
	if(is.null(sample.names)) sample.names <- removeExt(removeExt(removeExt(files)))
	if(readr) {
		OK <- requireNamespace("readr",quietly=TRUE)
		if(!OK) stop("readr package required but is not installed (or can't be loaded)")
	}

#	Read all files and store
	CountList <- list()
	ChrRleList <- list()
	LocusList <- list()
	ChrNames <- c()
	MaxLocus <- 1L
	for(i in 1:nsamples) {
		if(verbose) cat("Reading",files[i],"\n")
		if(readr)
			x <- as.data.frame(suppressWarnings(readr::read_tsv(files[i],col_names=FALSE,col_types="ci__ii",progress=FALSE)))
		else
			x <- read.delim(files[i],header=FALSE,colClasses=c("character","integer","NULL","NULL","integer","integer"))
		ChrRleList[[i]] <- rle(x[,1])
		LocusList[[i]] <- x[,2]
		CountList[[i]] <- as.matrix(x[,3:4])
		ChrNames <- unique(c(ChrNames,ChrRleList[[i]]$values))
	}

	if(verbose) cat("Hashing ...\n")

#	Convert rle values to integer
	for(i in 1:nsamples) ChrRleList[[i]]$values <- match(ChrRleList[[i]]$values,ChrNames)

#	Hash the genomic positions
	HashBase <- length(ChrNames)+1L
	HashList <- list()
	HashUnique <- c()
	for (i in 1:nsamples) {
		HashList[[i]] <- inverse.rle(ChrRleList[[i]]) / HashBase + LocusList[[i]]
		HashUnique <- unique(c(HashUnique,HashList[[i]]))
	}

	if(verbose) cat("Collating counts ...\n")

#	Merged count matrix
	counts <- matrix(0L,length(HashUnique),nsamples*2L)
	j <- 1:2
	for(i in 1:nsamples) {
		m <- match(HashList[[i]], HashUnique)
		counts[m,j] <- CountList[[i]]
		j <- j+2L
	}

#	Unhash
	Locus <- as.integer(HashUnique)
	Chr <- as.integer( (HashUnique-Locus) * HashBase + 0.5 )
	attr(Chr,"levels") <- ChrNames
	class(Chr) <- "factor"

#	Attach dimension names and form DGEList
	Sample2 <- rep(sample.names, each=2L)
	Methylation <- rep.int(c("Me","Un"), nsamples)
	colnames(counts) <- paste(Sample2, Methylation, sep="-")
	y <- DGEList(counts, genes=data.frame(Chr,Locus))
	row.names(y) <- paste(Chr,Locus,sep="-")
	y
}

