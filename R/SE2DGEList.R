SE2DGEList <- function(object)
#	Given any SummarizedExperiment data object, extract basic information needed
#	and convert it into a DGEList object
#	Yunshun Chen, Gordon Smyth
#	18 March 2020. Last modified 23 March 2020.
{
	if(!is(object,"SummarizedExperiment"))
		stop("object is not of the SummarizedExperiment class")

	if(!requireNamespace("SummarizedExperiment",quietly=TRUE))
		stop("SummarizedExperiment package required but is not installed (or can't be loaded)")

#	Check 'assays'
	if( !("counts" %in% SummarizedExperiment::assayNames(object)) ) stop("object doesn't contain counts")
	counts <- SummarizedExperiment::assay(object,"counts")

	if(!is.null(rownames(object))) rownames(counts) <- rownames(object)
	if(!is.null(colnames(object))) colnames(counts) <- colnames(object)

	genes <- samples <- NULL

#	Check 'colData'
	if(ncol(SummarizedExperiment::colData(object))){
		samples <- as.data.frame(SummarizedExperiment::colData(object))
	}

#	Check 'rowData'
	if(is(SummarizedExperiment::rowRanges(object), "GRanges")) 
		genes <- as.data.frame(SummarizedExperiment::rowRanges(object))
	else if(ncol(SummarizedExperiment::rowData(object)))
		genes <- as.data.frame(SummarizedExperiment::rowData(object))

	DGEList(counts=counts, samples=samples, genes=genes)
}
