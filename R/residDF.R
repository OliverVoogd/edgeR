.comboGroups <- function(truths) 
# Function that returns a list of vectors of indices,
# where each vector refers to the rows with the same
# combination of TRUE/FALSE values in 'truths'.
# 
# written by Aaron Lun
# Created 24 October 2014
{
#	Integer packing will only work for 31 libraries at a time.	
	assembly <- list()
	collected <- 0L
	step <- 31L
	bits <- as.integer(2^(1:step-1L))

	while (collected < ncol(truths)) {
		upper <- pmin(ncol(truths) - collected, step)
		keys <- t(truths[,collected+1:upper,drop=FALSE]) * bits[1:upper]
		assembly[[length(assembly)+1L]] <- as.integer(colSums(keys))
		collected <- collected + step
	}

#	Figuring out the unique components.
	o <- do.call(order, assembly)
	nr <- nrow(truths)
	is.different <- logical(nr)
	for (i in 1:length(assembly)) { 
		is.different <- is.different | c(TRUE, diff(assembly[[i]][o])!=0L)
	}
	first.of.each <- which(is.different)
	last.of.each <- c(first.of.each[-1]-1L, nr)

#	Returning the groups.
	output <- list()
	for (u in 1:length(first.of.each)) {
		output[[u]] <- o[first.of.each[u]:last.of.each[u]]
	}
	return(output)
}

.residDF <- function(zero, design)
#	Effective residual degrees of freedom after adjusting for exact zeros
#	Gordon Smyth and Aaron Lun
#	Created 6 Jan 2014.  Last modified 2 Jun 2020
{
	nlib <- ncol(zero)
	ncoef <- ncol(design)
	nzero <- as.integer(rowSums(zero))

#	Default is no zero
	DF <- rep_len(nlib-ncoef,length(nzero))

#	All zero case
	DF[nzero==nlib] <- 0L

#	Anything in between?
	somezero <- nzero>0L & nzero<nlib
	if(any(somezero)) {
		zero2 <- zero[somezero,,drop=FALSE]
		groupings <- .comboGroups(zero2)

#		Identifying the true residual d.f. for each of these rows.			
		DF2 <- nlib-nzero[somezero]
		for (u in 1:length(groupings)) {
			i <- groupings[[u]]
			zeroi <- zero2[i[1],]
			DF2[i] <- DF2[i]-qr(design[!zeroi,,drop=FALSE])$rank
		}
		DF2 <- pmax(DF2, 0L)
		DF[somezero] <- DF2
	}

	DF
}
