readDGE <- function(files,path=NULL,columns=c(1,2),group=NULL,labels=NULL,...) 
#	Read and collate a set of count data files, each file containing counts for one library
#	Created 2006.  Last modified 8 Sep 2021.
{
#	Create data.frame to hold sample information
	x <- list()
	if(is.data.frame(files)) {
		x$samples <- files
		if(!is.null(labels)) row.names(x$samples) <- labels
		if(is.null(x$samples$files)) stop("file names not found")
		x$samples$files <- as.character(x$samples$files)
	} else {
		files <- as.character(files)
		if(is.null(labels)) labels <- removeExt(files)
		x$samples <- data.frame(files=files,row.names=labels,stringsAsFactors=FALSE)
	}
	nfiles <- nrow(x$samples)

#	Set group factor
	if(!is.null(group)) x$samples$group <- group
	if(is.null(x$samples$group)) x$samples$group <- rep_len(1L,nfiles)
	x$samples$group <- as.factor(x$samples$group)

#	Read files into a list
	d <- taglist <- list()
	for (fn in x$samples$files) {
		if(!is.null(path)) fn <- file.path(path,fn)
		d[[fn]] <- read.delim(fn,...,stringsAsFactors=FALSE)
		taglist[[fn]] <- as.character(d[[fn]][,columns[1]])
		if(anyDuplicated(taglist[[fn]])) {
			ndup <- sum(duplicated(taglist[[fn]]))
			stop("There are ",ndup," repeated row names in ",fn,". Row names must be unique.")
		}
	}

#	Collate counts for unique tags
	tags <- unique(unlist(taglist))
	ntags <- length(tags)
	x$counts <- matrix(0,ntags,nfiles)
	dimnames(x$counts) <- list(Tags=tags,Samples=row.names(x$samples))
	for (i in 1:nfiles) {
		aa <- match(taglist[[i]],tags)
		x$counts[aa,i] <- d[[i]][,columns[2]]
	}

#	Alert user if htseq-style meta genes found
	MetaTags <- grep("^_",tags,value=TRUE)
	if(length(MetaTags)) message("Meta tags detected: ",paste(MetaTags,collapse=", "))

#	Enter library sizes and norm factors
	x$samples$lib.size <- colSums(x$counts)
	x$samples$norm.factors <- 1

	new("DGEList",x)
}
