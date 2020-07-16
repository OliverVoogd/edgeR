# S3 as.matrix method

as.matrix.DGEList <- function(x,...) as.matrix(x$counts)

# S3 as.data.frame method

as.data.frame.DGEExact <- as.data.frame.DGELRT <- function(x,row.names=NULL,...)
{
	if(is.null(x$genes)) {
		if(!is.null(row.names)) row.names(x$table) <- row.names
		x$table		
	} else
		data.frame(x$genes,x$table,row.names=row.names,check.rows=FALSE,check.names=FALSE,stringsAsFactors=FALSE)
}
as.data.frame.TopTags <- function(x,row.names=NULL,...)
{
	if(!is.null(row.names)) row.names(x$table) <- row.names
	x$table
}

# S3 dim methods
# These enable nrow() and ncol() as well

dim.DGEList <- function(x) if(is.null(x$counts)) c(0,0) else dim(as.matrix(x$counts))
dim.DGEGLM <- function(x) if(is.null(x$coefficients)) c(0,0) else dim(as.matrix(x$coefficients))
dim.DGEExact <- dim.TopTags <- dim.DGELRT <- function(x) if(is.null(x$table)) c(0,0) else dim(as.matrix(x$table))

# S3 length methods
# These methods have been removed so that edgeR objects will return list length instead of matrix length

#length.DGEList <- length.DGEExact <- length.TopTags <- length.DGEGLM <- length.DGELRT <- function(x) prod(dim(x))

# S3 dimnames methods
# These enable rownames() and colnames() as well

dimnames.DGEList <- function(x) dimnames(x$counts)
dimnames.DGEGLM <- function(x) dimnames(x$coefficients)
dimnames.DGEExact <- dimnames.DGELRT <- dimnames.TopTags <- function(x) dimnames(x$table)

# S3 dimnames<- methods
# These enable rownames()<- and colnames()<- as well

assign("dimnames<-.DGEList",function(x,value)
{
	dimnames(x$counts) <- value
	if(!is.null(x$samples)) row.names(x$samples) <- value[[2]]
	if(!is.null(x$genes)) row.names(x$genes) <- value[[1]]
	x
})

assign("dimnames<-.DGEExact",function(x,value)
{
	dimnames(x$table) <- value
	if(!is.null(x$genes)) row.names(x$genes) <- value[[1]]
	x
})

assign("dimnames<-.DGEGLM",function(x,value)
{
	dimnames(x$coefficients) <- value
	if(!is.null(x$unshrunk.coefficients)) dimnames(x$unshrunk.coefficients) <- value
	if(!is.null(x$fitted.values)) rownames(x$fitted.values) <- value[[1]]
	if(!is.null(x$counts)) rownames(x$fitted.values) <- value[[1]]
	if(!is.null(x$genes)) row.names(x$genes) <- value[[1]]
	x
})

assign("dimnames<-.DGELRT",function(x,value)
#	4 June 2015
{
	dimnames(x$table) <- value
	if(!is.null(x$coefficients)) rownames(x$coefficients) <- value[[1]]
	if(!is.null(x$unshrunk.coefficients)) rownames(x$unshrunk.coefficients) <- value[[1]]
	if(!is.null(x$fitted.values)) rownames(x$fitted.values) <- value[[1]]
	if(!is.null(x$counts)) rownames(x$fitted.values) <- value[[1]]
	if(!is.null(x$genes)) row.names(x$genes) <- value[[1]]
	x
})

# S3 head and tail methods

head.DGEList <- head.DGEExact <- head.DGEGLM <- head.DGELRT <- head.TopTags <-
function (x, n = 6L, ...) 
{
	stopifnot(length(n) == 1L)
	n <- if (n < 0L) 
		max(nrow(x) + n, 0L)
	else
		min(n, nrow(x))
	x[seq_len(n),]
}

tail.DGEList <- tail.DGEExact <- tail.DGEGLM <- tail.DGELRT <- tail.TopTags <-
function (x, n = 6L, ...) 
{
	stopifnot(length(n) == 1L)
	nrx <- nrow(x)
	n <- if (n < 0L) 
		max(nrx + n, 0L)
	else
		min(n, nrx)
	x[seq.int(to = nrx, length.out = n),]
}
