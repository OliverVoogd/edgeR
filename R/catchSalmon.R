catchSalmon <- function(paths,verbose=TRUE)
#	Read transcriptwise counts and bootstrap samples from Salmon output
#	Use bootstrap samples to estimate overdispersion of transcriptwise counts
#	Gordon Smyth
#	Created 1 April 2018. Last modified 28 Aug 2019.
{
	NSamples <- length(paths)

#	Use jsonlite and readr packages for reading
	OK <- requireNamespace("jsonlite",quietly=TRUE)
	if(!OK) stop("jsonlite package required but is not installed (or can't be loaded)")
	OK <- requireNamespace("readr",quietly=TRUE)
	if(!OK) stop("readr package required but is not installed (or can't be loaded)")

#	Accumulate counts and CV^2 of bootstrap counts for each sample
	for (j in 1L:NSamples) {
		if(verbose) cat("Reading ",paths[j],", ",sep="")

#		File locations
		MetaFile <- file.path(paths[j],"aux_info","meta_info.json")
		QuantFile <- file.path(paths[j],"quant.sf")
		BootFile <- file.path(paths[j],"aux_info","bootstrap","bootstraps.gz")
		if(!file.exists(QuantFile)) stop("quant.sf file not found at specified path")

#		Meta information
		Meta <- jsonlite::fromJSON(MetaFile)
		NTx <- Meta$num_targets
		if(is.null(NTx)) NTx <- Meta$num_valid_targets
		if(is.null(NTx)) stop("Can't find number of targets")
		NBoot <- Meta$num_bootstraps
		if(is.null(NBoot)) stop("Can't find number of bootstraps")
		if(verbose) cat(NTx,"transcripts,",NBoot,"bootstraps\n")

#		Read counts
		if(j == 1L) {
			Counts <- matrix(0,NTx,NSamples)
			DF <- rep_len(0L,NTx)
			OverDisp <- rep_len(0,NTx)
			Quant1 <- suppressWarnings(readr::read_tsv(QuantFile,col_types="cdd_d",progress=FALSE))
			Counts[,1L] <- Quant1$NumReads	
		} else {
			Quant <- suppressWarnings(readr::read_tsv(QuantFile,col_types="____d",progress=FALSE))
			Counts[,j] <- Quant$NumReads
		}

#		Bootstrap samples
		if(NBoot > 0L) {
			BootFileCon <- gzcon(file(BootFile,open="rb"))
			Boot <- readBin(BootFileCon,what="double",n=NTx*NBoot)
			close(BootFileCon)
			dim(Boot) <- c(NTx,NBoot)
			M <- rowMeans(Boot)
			i <- (M > 0)
			OverDisp[i] <- OverDisp[i] + rowSums((Boot[i,]-M[i])^2) / M[i]
			DF[i] <- DF[i]+NBoot-1L
		}
	}

#	Estimate overdispersion for each transcript
	i <- (DF > 0L)
	if(sum(i) > 0L) {
		OverDisp[i] <- OverDisp[i] / DF[i]
#		Apply a limited amount of moderation
		DFMedian <- median(DF[i])
		DFPrior <- 3
		OverDispPrior <- median(OverDisp[i]) / qf(0.5,df1=DFMedian,df2=DFPrior)
		if(OverDispPrior < 1) OverDispPrior <- 1
		OverDisp[i] <- (DFPrior * OverDispPrior + DF[i]*OverDisp[i]) / (DFPrior + DF[i])
		OverDisp <- pmax(OverDisp,1)
		OverDisp[!i] <- OverDispPrior
	} else {
		OverDisp[] <- NA_real_
		OverDispPrior <- NA_real_
	}

#	Prepare output
	Quant1 <- as.data.frame(Quant1,stringsAsFactors=FALSE)
	dimnames(Counts) <- list(Quant1$Name,paths)
	row.names(Quant1) <- Quant1$Name
	Quant1$Name <- NULL
	Quant1$TPM <- Quant1$NumReads <- NULL
	Quant1$Overdispersion <- OverDisp

	list(counts=Counts,annotation=Quant1,overdispersion.prior=OverDispPrior)
}
