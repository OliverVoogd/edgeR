# Split the data according to group

splitIntoGroups <- function(y, ...)
UseMethod("splitIntoGroups")

splitIntoGroups.DGEList <- function(y, ...)
# Yunshun Chen. Created 18 March 2016.
{
	group <- y$samples$group
	splitIntoGroups(y$counts, group=group)
}

splitIntoGroups.default <- function(y, group=NULL, ...) 
# Written by Davis McCarthy, February 2009, idea suggested by Mark Robinson
# Last modified 2 Jun 2020.
{
#	Check y
	y <- as.matrix(y)
	ntags <- nrow(y)
	nlibs <- ncol(y)

#	Check group
	if(is.null(group)) group <- rep_len(1L, nlibs)
	if(length(group)!=nlibs) stop("Incorrect length of group.")
	group <- dropEmptyLevels(group)
	
	out <- lapply(split(t(y), group), FUN=function(u) matrix(u, nrow=ntags, byrow=TRUE))
	out
}