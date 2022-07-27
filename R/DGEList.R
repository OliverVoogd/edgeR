DGEList <- function(counts=matrix(0,0,0), lib.size=colSums(counts), norm.factors=rep(1,ncol(counts)), samples=NULL, group=NULL, genes=NULL, remove.zeros=FALSE) 
#	Construct DGEList object from components, with some checking
#	Created 28 Sep 2008. Last modified 16 Nov 2020.
{
#	Check whether counts is a data.frame
	if(is.data.frame(counts)) {
		cl <- vapply(counts,class,FUN.VALUE="")
		IsNumeric <- (cl %in% c("integer","numeric"))
		a <- which(!IsNumeric)
		if(length(a)) {
			if(length(a)==1L && a[1]==1L) {
				stop(
"The count matrix is a data.frame instead of a matrix and the first column is of class ",cl[1],"\n",
"  instead of being numeric. Was the first column intended to contain geneids?"
				)
			} else {
				stop(
"The count matrix is a data.frame instead of a matrix and ",length(a)," columns are non-numeric.\n",
"  Should these columns be gene annotation instead of counts?"
				)
			}
		}
	}

#	Check counts
	counts <- as.matrix(counts)
	if( !any(typeof(counts) == c("integer","double")) ) stop("non-numeric values found in counts")
	nlib <- ncol(counts)
	ntags <- nrow(counts)
	if(nlib>0L && is.null(colnames(counts))) colnames(counts) <- paste0("Sample",1L:nlib)
	if(ntags>0L && is.null(rownames(counts))) rownames(counts) <- 1L:ntags
	.isAllZero(counts) # don't really care about all-zeroes, but do want to protect against NA's, negative values.

#	Check lib.size
	if(is.null(lib.size)) {
		lib.size <- colSums(counts)
		if(min(lib.size) <= 0) warning("library size of zero detected")
	} else {
		if(!is.numeric(lib.size)) stop("'lib.size' must be numeric")
		if(nlib != length(lib.size)) stop("length of 'lib.size' must equal number of columns in 'counts'")
		minlibsize <- min(lib.size)
		if(is.na(minlibsize)) stop("NA library sizes not allowed")
		if(minlibsize < 0) stop("negative library sizes not permitted")
		if(minlibsize == 0) {
			if(any(lib.size==0 & colSums(counts)>0)) stop("library size set to zero but counts are nonzero")
			warning("library size of zero detected")
		}
	}

#	Check norm.factors
	if(is.null(norm.factors)) {
		norm.factors <- rep_len(1,ncol(counts))
	} else {
		if(!is.numeric(norm.factors)) stop("'lib.size' must be numeric")
		if(!identical(nlib,length(norm.factors))) stop("Length of 'norm.factors' must equal number of columns in 'counts'")
		minnf <- min(norm.factors)
		if(is.na(minnf)) stop("NA norm factors not allowed")
		if(minnf <= 0) stop("norm factors should be positive")
		if( abs(prod(norm.factors) - 1) > 1e-6 ) warning("norm factors don't multiply to 1")
	}

#	Check samples
	if(!is.null(samples)) {
		samples <- as.data.frame(samples)
		if(nlib != nrow(samples)) stop("Number of rows in 'samples' must equal number of columns in 'counts'")
	}

#	Get group from samples if appropriate
	if(is.null(group) && !is.null(samples$group)) {
		group <- samples$group
		samples$group <- NULL
	}

#	Check group
	if(is.null(group)) {
		group <- rep_len(1L,nlib)
		levels(group) <- "1"
		class(group) <- "factor"
	} else {
		if(length(group) != nlib) stop("Length of 'group' must equal number of columns in 'counts'")
		group <- dropEmptyLevels(group)
	}

#	Make data frame of sample information
	sam <- data.frame(group=group,lib.size=lib.size,norm.factors=norm.factors)
	if(!is.null(samples)) sam <- data.frame(sam, samples)
	samples <- sam
	if(anyDuplicated(colnames(counts))) {
		message("Repeated column names found in count matrix")
		row.names(samples) <- 1L:nlib
	} else 
		row.names(samples) <- colnames(counts)

#	Make object
	x <- new("DGEList",list(counts=counts,samples=samples))

#	Add data frame of gene information
	if(!is.null(genes)) {
		genes <- as.data.frame(genes, stringsAsFactors=FALSE)
		if(nrow(genes) != ntags) stop("Counts and genes have different numbers of rows")
		if(anyDuplicated(row.names(counts)))
			warning("Count matrix has duplicated rownames",call.=FALSE)
		else 
			row.names(genes) <- row.names(counts)
		x$genes <- genes
	}

#	Remove rows with all zeros
	if(remove.zeros) {
		all.zeros <- rowSums(counts>0,na.rm=TRUE)==0
		if(any(all.zeros)) {
			x <- x[!all.zeros,]
			message("Removing ",sum(all.zeros)," rows with all zero counts")
		}
	}

#	x$offset <- expandAsMatrix(getOffset(x),dim(counts))
#	x$weights <- matrix(1,ntags,nlib)

	x
}

.isAllZero <- function(y) 
# Check whether all counts are zero.
# Also checks and stops with an informative error message if negative, NA or infinite counts are present.
{
	if (!length(y)) return(FALSE)
	check.range <- range(y)
	if (is.na(check.range[1])) stop("NA counts not allowed", call.=FALSE)
	if (check.range[1] < 0) stop("Negative counts not allowed", call.=FALSE)
	if (is.infinite(check.range[2])) stop("Infinite counts not allowed", call.=FALSE)
	check.range[2]==0
}
