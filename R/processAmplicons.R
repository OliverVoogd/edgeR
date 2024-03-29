processAmplicons = function(readfile, readfile2=NULL, barcodefile, hairpinfile,
                    allowMismatch=FALSE, barcodeMismatchBase = 1, hairpinMismatchBase = 2,
                    dualIndexForwardRead=FALSE, verbose = FALSE, barcodesInHeader = FALSE,
					hairpinBeforeBarcode=FALSE,
                    plotPositions=FALSE) {

  ## a simple check for the existance of the given files
  checkFileExistence = function(readfilenames){
    for(i in 1:length(readfilenames)){
      if (!file.exists(readfilenames[i]))
        stop(paste("Read file ", readfilenames[i], " doesn't exist. \n", sep=""))
    }
  }

  # Check file existence
  numfiles = length(readfile)
  checkFileExistence(readfile);
  IsPairedReads = FALSE;
  if (!is.null(readfile2)) {
    IsPairedReads = TRUE;
    if (numfiles != length(readfile2))
      stop("readfile and readfile2 length should match each other.")
    checkFileExistence(readfile2);
  }

  IsDualIndexingOnForwardRead = dualIndexForwardRead

  # check for existance of barcode and hairpin file
  if (!file.exists(barcodefile))
    stop("Barcode file doesn't exist.\n")
  if (!file.exists(hairpinfile))
    stop("Hairpin file doesn't exist.\n")

  # Find the length of the first read
  reads <- file(readfile[1], "rt");
  first_read <- readLines(reads, 2)
  readlength <- nchar(first_read[2])
  close(reads)

  # Validate barcodes
  barcodes <- read.table(barcodefile, header=TRUE, sep="\t")
  barcodeIDIndex = which(colnames(barcodes) == 'ID')

  if (length(barcodeIDIndex) < 1)
    stop("Can't find column ID in ", barcodefile)

  barcodeseqIndex = which(colnames(barcodes) == 'Sequences')
  if (length(barcodeseqIndex) < 1)
    stop("Can't find column Sequences in ", barcodefile)

  barcodeIDs <- as.character(barcodes[, barcodeIDIndex])
  barcodeseqs <- as.character(barcodes[, barcodeseqIndex])
  barcodelength <- nchar(barcodeseqs[[1]])
  barcode2length <- 0
  barcodelengthReverse <- 0

  if (anyDuplicated(barcodeIDs))
    stop("There are duplicate barcode IDs.\n")

  if ((min(nchar(barcodeseqs)) != barcodelength) || (max(nchar(barcodeseqs)) != barcodelength)) ## barcodes all must be same length
    stop(paste("Barcode sequence length is set to ", barcodelength, ", there are barcode sequence not with specified length.\n", sep=""))

  if (IsPairedReads) {
    barcodeseqRevIndex = which(colnames(barcodes) == 'SequencesReverse') ## finds the index of reverse barcodes
    if (length(barcodeseqRevIndex) < 1)
      stop("Can't find column SequencesReverse in ", barcodefile)
    barcodeseqsReverse <- as.character(barcodes[, barcodeseqRevIndex]) ## finds the reverse barcodes
    barcodelengthReverse <- nchar(barcodeseqsReverse[[1]])
    if ((min(nchar(barcodeseqsReverse)) != barcodelengthReverse) || (max(nchar(barcodeseqsReverse)) != barcodelengthReverse))
      stop(paste("Reverse barcode sequence length is set to ", barcodelength, ", there are reverse barcode sequence not in specified length.\n", sep=""))

    concatenatedBarcodeseqs = paste(barcodeseqs, barcodeseqsReverse, sep="")
    if (anyDuplicated(concatenatedBarcodeseqs))
      stop("There are duplicate forward/reverse barcode sequences.\n")

  } else if (IsDualIndexingOnForwardRead) {
    barcodeseq2Index = which(colnames(barcodes) == 'Sequences2')
    if (length(barcodeseq2Index) < 1)
      stop("Can't find column Sequences2 in ", barcodefile)

    barcode2seqs <- as.character(barcodes[, barcodeseq2Index])
    barcode2length <- nchar(barcode2seqs[[1]])

    if ((min(nchar(barcode2seqs)) != barcode2length) || (max(nchar(barcode2seqs)) != barcode2length))
      stop(paste("Forward barcode2 sequence length is set to ", barcode2length, ", there are barcode2 sequence not in specified length.\n", sep=""))
    concatenatedBarcodeseqs = paste(barcodeseqs, barcode2seqs, sep="")
    if (anyDuplicated(concatenatedBarcodeseqs))
      stop("There are duplicate barcode/barcode2 sequences.\n")

  } else {
    if (anyDuplicated(barcodeseqs)) stop("There are duplicate barcode sequences.\n")
  }

  # Validate hairpins
  hairpins <- read.table(hairpinfile, header=TRUE, sep="\t"); ## read in the hairpin table

  hairpinIDIndex = which(colnames(hairpins) == 'ID')
  if (length(hairpinIDIndex) < 1)
    stop("Can't find column ID in ", hairpinfile)
  hairpinIDs <- as.character(hairpins[, hairpinIDIndex])

  hairpinSeqIndex = which(colnames(hairpins) == 'Sequences')
  if (length(hairpinSeqIndex) < 1)
    stop("Can't find column Sequences in ", hairpinfile)
  hairpinseqs <- as.character(hairpins[, hairpinSeqIndex])
  hairpinlength <- nchar(hairpinseqs[[1]])

  if ((min(nchar(hairpinseqs)) != hairpinlength) || (max(nchar(hairpinseqs)) != hairpinlength))
    stop(paste("Hairpin sequence length is set to ", hairpinlength, ", there are hairpin sequences not with specified length.\n", sep=""))
  if (anyDuplicated(hairpinseqs))
    stop("There are duplicate hairpin sequences.\n")
  if (anyDuplicated(hairpinIDs))
    stop("There are duplicate hairpin IDs.\n")

  # validate mismatch/shifting input parameters
  if (allowMismatch) {
    if ((barcodeMismatchBase < 0) || (barcodeMismatchBase > 2))
      stop("To allow mismatch in barcode sequence, please input a non-negative barcodeMismatchBase no greater than than 2. ")
    if ((hairpinMismatchBase < 0) || (hairpinMismatchBase > 4))
      stop("To allow mismatch in hairpin sequence, please input a non-negative hairpinMismatchBase no greater than than 4. ")
  }

  tempbarcodefile <- paste("Barcode", as.character(Sys.getpid()), "temp.txt", sep = "_")
  on.exit({ if (file.exists(tempbarcodefile)) { file.remove(tempbarcodefile) }}, add=TRUE)

  if (IsPairedReads) {
    bothBarcodeSeqs = cbind(barcodeseqs, barcodeseqsReverse)
    write.table(bothBarcodeSeqs, file=tempbarcodefile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE);
  } else if (IsDualIndexingOnForwardRead) {
    bothBarcodeSeqs = cbind(barcodeseqs, barcode2seqs)
    write.table(bothBarcodeSeqs, file=tempbarcodefile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE);
  } else {
    write.table(barcodeseqs, file=tempbarcodefile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE);
  }

  temphairpinfile <- paste("Hairpin", as.character(Sys.getpid()), "temp.txt", sep = "_")
  on.exit({ if (file.exists(temphairpinfile)) { file.remove(temphairpinfile) }}, add=TRUE) ## remove file on exit
  write.table(hairpinseqs, file=temphairpinfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE);

  ## generate a unique output file for C output storage
  tempoutfile <- paste("ReadcountSummary", as.character(Sys.getpid()), "output.txt", sep = "_")
  on.exit({ if (file.exists(tempoutfile)) { file.remove(tempoutfile) }}, add=TRUE)

  ## generate unique files for graphing the distribution of positive barcode and hairpin read locations
  tempbarcodeposfile <- paste("BarcodePosition", as.character(Sys.getpid()), "Summary.txt", sep="_")
  on.exit({if (file.exists(tempbarcodeposfile)) { file.remove(tempbarcodeposfile) }}, add=TRUE)

  tempbarcode2posfile <- paste("Barcode2Position", as.character(Sys.getpid()), "Summary.txt", sep="_")
  on.exit({if (file.exists(tempbarcode2posfile)) { file.remove(tempbarcode2posfile) }}, add=TRUE)

  temphairpinposfile <- paste("HairpinPosition", as.character(Sys.getpid()), "Summary.txt", sep="_")
  on.exit({if (file.exists(temphairpinposfile)) {file.remove(temphairpinposfile) }}, add=TRUE)

  tryCatch({
    if (!IsPairedReads) {
       readfile2 = rep("DummyReadfile.fastq", numfiles) ## create fake files if only one file given
       barcodelengthReverse = 0
    }
    if (!IsDualIndexingOnForwardRead) {
       barcode2length = 0
    }

    ## calls the c function .cxx_processHairpinReads which fills data into tempoutfile
    .C(.cxx_processHairpinReads, as.integer(IsPairedReads), as.integer(IsDualIndexingOnForwardRead),
	     as.character(readfile), as.character(readfile2), as.integer(numfiles),
       as.character(tempbarcodefile), as.character(temphairpinfile),
	     as.integer(barcodelength), as.integer(barcode2length), as.integer(barcodelengthReverse),
	     as.integer(hairpinlength),
       as.integer(allowMismatch), as.integer(barcodeMismatchBase), as.integer(hairpinMismatchBase),
       as.character(tempoutfile), as.integer(verbose), as.integer(barcodesInHeader), as.integer(plotPositions),
	     as.character(tempbarcodeposfile), as.character(tempbarcode2posfile), as.character(temphairpinposfile),
		 as.integer(hairpinBeforeBarcode))

    ## retrive all of the calculated data
    hairpinReadsSummary <- read.table(tempoutfile, sep="\t", header=FALSE)

    ## if plotPositions is true, plot all of the relevant data
    if (plotPositions) {
      numberOfPlots <- 1

      oneRowFrame <- data.frame(read_position = c(0), counts = c(0))
      barcodePositionSummary <- read.table(tempbarcodeposfile, sep="\t", header=FALSE)
      barcodePositionSummary <- data.frame(read_position= c(1:nrow(barcodePositionSummary)),
                                           counts = barcodePositionSummary$V1)
      barcodePositionSummary <- rbind(oneRowFrame, barcodePositionSummary)

      if (IsPairedReads || IsDualIndexingOnForwardRead) {
        numberOfPlots <- 2

        barcode2PositionSummary <- read.table(tempbarcode2posfile, sep="\t", header=FALSE)
        barcode2PositionSummary <- data.frame(read_position=c(1:nrow(barcode2PositionSummary)),
                                              counts=barcode2PositionSummary$V1)
        barcode2PositionSummary <- rbind(oneRowFrame, barcode2PositionSummary)
      }

      hairpinPositionSummary <- read.table(temphairpinposfile, sep="\t", header=FALSE)
      hairpinPositionSummary <- data.frame(read_position = c(1:nrow(hairpinPositionSummary)),
                                           counts = hairpinPositionSummary$V1)
      hairpinPositionSummary <- rbind(oneRowFrame, hairpinPositionSummary)

      par(mfrow=c(1, numberOfPlots))
      plotColours <- c("firebrick", "navyblue")
      legendNames <- c("Barcodes", "Hairpins")
      plot(NULL,
           xlim=c(0, max(max(barcodePositionSummary$read_position), max(hairpinPositionSummary$read_position))),
           ylim=c(0, 1.2 * max(max(barcodePositionSummary$counts, max(hairpinPositionSummary$counts)))),
           ylab = "Sequence Counts",
           xlab = "Read Position",
           main = "Barcode & Hairpin Position")
      polygon(barcodePositionSummary$read_position, barcodePositionSummary$counts, col="firebrick", border="firebrick")
      polygon(hairpinPositionSummary$read_position, hairpinPositionSummary$counts, col="navyblue", border="navyblue")

      legend(x=max(barcodePositionSummary$read_position),
             y= 1.2 * max(barcodePositionSummary$counts),
             legend=legendNames,
             col=plotColours,
             fill=plotColours,
             xjust=1)


      if (IsPairedReads || IsDualIndexingOnForwardRead) {
        title <- "Paired Read Barcode Positions"
        if (IsDualIndexingOnForwardRead) {
         title <- "Dual Indexed Barcode Positions"
        }

        plot(NULL,
             xlim=c(0, max(barcode2PositionSummary$read_position)),
             ylim=c(0, 1.2 * max(barcode2PositionSummary$counts)),
             ylab = "Sequence Counts",
             xlab = "Read Position",
             main = title)
        polygon(barcode2PositionSummary$read_position, barcode2PositionSummary$counts, col="firebrick", border="firebrick")
      }
    }
  }, error = function(err) {print(paste("ERROR MESSAGE:  ",err))}
  )


  if (exists("hairpinReadsSummary")) {
    if (nrow(hairpinReadsSummary) != length(hairpinIDs))
      stop("Number of hairpins from result count matrix doesn't agree with given hairpin list. ")
    if (ncol(hairpinReadsSummary) != length(barcodeIDs))
      stop("Number of barcodes from result count matrix doesn't agree with given barcode list. ")
    colnames(hairpinReadsSummary) = barcodeIDs
    rownames(hairpinReadsSummary) = hairpinIDs

    x <- edgeR::DGEList(counts = hairpinReadsSummary, genes = hairpins)
    if(!is.null(barcodes$group)) {
      x$samples = cbind("ID"=barcodes$ID, "lib.size"=x$samples$lib.size,
                       "norm.factors"=x$samples$norm.factors,
                       barcodes[,-match(c("ID","Sequences"), colnames(barcodes))])
    } else {
      x$samples = cbind("ID"=barcodes$ID, x$samples, barcodes[,-match(c("ID","Sequences"), colnames(barcodes))])
    }
  } else {
    stop("An error occured in processHairpinReads.")
  }

  # remove all of these??
  #if (file.exists(tempbarcodefile)) { file.remove(tempbarcodefile) }
  #if (file.exists(temphairpinfile)) { file.remove(temphairpinfile) }
  #if (file.exists(tempoutfile)) { file.remove(tempoutfile) }
  #if (file.exists(tempbarcodeposfile)) { file.remove(tempbarcodeposfile) }
  #if (file.exists(tempbarcode2posfile)) { file.remove(tempbarcode2posfile) }
  #if (file.exists(temphairpinposfile)) {file.remove(temphairpinposfile) }
  return(x)
}
