read10X <- function(mtx=NULL,genes=NULL,barcodes=NULL,path=".",DGEList=TRUE)
#	Read 10X Genomics Matrix Exchange Format files created by CellRanger
#	Gordon Smyth
#	Created 10 Jan 2018. Last modified 24 Apr 2019.
{
#	Get file names
	if(is.null(mtx) || is.null(genes) || is.null(barcodes)) {
		files <- dir(path)
		if(is.null(mtx)) {
			if("matrix.mtx.gz" %in% files) mtx <- "matrix.mtx.gz"
			if("matrix.mtx" %in% files) mtx <- "matrix.mtx"
			if(is.null(mtx)) stop("Can't find matrix.mtx file, please specify filename explicitly")
		}
		if(is.null(genes)) {
			if("genes.tsv.gz" %in% files) genes <- "genes.tsv.gz"
			if("genes.tsv" %in% files) genes <- "genes.tsv"
			if("features.tsv.gz" %in% files) genes <- "features.tsv.gz"
			if("features.tsv" %in% files) genes <- "features.tsv"
			if(is.null(genes)) stop("Can't find genes.tsv or features.tsv, please specify filename explicitly")
		}
		if(is.null(barcodes)) {
			if("barcodes.tsv" %in% files) barcodes <- "barcodes.tsv"
			if("barcodes.tsv.gz" %in% files) barcodes <- "barcodes.tsv.gz"
		}
	}

#	Add path
	mtx <- file.path(path,mtx)
	genes <- file.path(path,genes)
	if(!is.null(barcodes)) barcodes <- file.path(path,barcodes)

#	Fetch header info for checking
	N <- scan(mtx,skip=2,what=0L,sep=" ",nmax=3,quiet=TRUE)
	ngenes <- N[1]
	ncells <- N[2]
	nmtx <- N[3]

#	Read gene Ids
	Genes <- read.table(genes,header=FALSE,comment.char="",sep="\t",row.names=1,colClasses="character")
	if(nrow(Genes) != ngenes) stop("Number of feature IDs doesn't agree with header information in mtx file")
	names(Genes)[1] <- "Symbol"
	if(ncol(Genes) > 1L) names(Genes)[2] <- "Type"

#	Read mtx file of counts
	m <- read.table(mtx,skip=3,header=FALSE,comment.char="",sep=" ",colClasses="integer",nrows=nmtx)

#	Convert Market Exchange Format to ordinary matrix
	y <- matrix(0L,ngenes,ncells)
	i <- m[,1]+(m[,2]-1L)*ngenes
	y[i] <- m[,3]
	dimnames(y) <- list(Gene=row.names(Genes),Cell=1:ncells)

#	Optionally read barcodes
	if(is.null(barcodes)) {
		Samples <- NULL
	} else {
		Barcodes <- scan(barcodes,what="",quiet=TRUE)
		if(length(Barcodes) != ncells) stop("Number of barcodes doesn't agree with header information in mtx file")
		Samples <- data.frame(Barcode=Barcodes)
	}

	if(DGEList) {
		DGEList(count=y,genes=Genes,samples=Samples)
	} else {
		list(counts=y,samples=Samples,genes=Genes)
	}
}
