exactTest <- function(object, pair=1:2, dispersion="auto", rejection.region="doubletail", big.count=900, prior.count=0.125)
#	Calculates exact p-values for the differential expression levels of tags in the two groups being compared.
#	Davis McCarthy, Gordon Smyth.
#	Created September 2009. Last modified 8 July 2012.
{
#	Check input
	if(!is(object,"DGEList")) stop("Currently only supports DGEList objects as the object argument.")
	if(length(pair)!=2) stop("Pair must be of length 2.")
	rejection.region <- match.arg(rejection.region,c("doubletail","deviance","smallp"))

#	Get group names
	group <- as.factor(object$samples$group)
	levs.group <- levels(group)
	if(is.numeric(pair))
		pair <- levs.group[pair]
	else
		pair <- as.character(pair)	
	if(!all(pair %in% levs.group)) stop("At least one element of given pair is not a group.\n Groups are: ", paste(levs.group, collapse=" "))

#	Get dispersion vector
	if(is.null(dispersion)) dispersion <- "auto"
	if(is.character(dispersion)) {
		dispersion <- match.arg(dispersion,c("auto","common","trended","tagwise"))
		dispersion <- switch(dispersion,
			"common"=object$common.dispersion,
			"trended"=object$trended.dispersion,
			"tagwise"=object$tagwise.dispersion,
			"auto"=getDispersion(object)
		)
		if(is.null(dispersion)) stop("specified dispersion not found in object")
		if(is.na(dispersion[1])) stop("dispersion is NA")
	}
	ldisp <- length(dispersion)
	ntags <- nrow(object$counts)
	if(ldisp!=1 && ldisp!=ntags) stop("Dispersion provided by user must have length either 1 or the number of tags in the DGEList object.")
	if(ldisp==1) dispersion <- rep(dispersion,ntags)

#	Reduce to two groups
	group <- as.character(group)
	j <- group %in% pair
	y <- object$counts[,j,drop=FALSE]
	lib.size <- object$samples$lib.size[j]
	norm.factors <- object$samples$norm.factors[j]
	group <- group[j]
	if(is.null(rownames(y))) rownames(y) <- paste("tag",1:ntags,sep=".")

#	Normalized library sizes
	lib.size <- lib.size * norm.factors
	offset <- log(lib.size)
	lib.size.average <- exp(mean(offset))

#	logFC
	prior.count <- prior.count*lib.size/mean(lib.size)
	offset.aug <- log(lib.size+2*prior.count)
	j1 <- group==pair[1]
	n1 <- sum(j1)
	if(n1==0) stop("No libraries for",pair[1])
	y1 <- y[,j1,drop=FALSE]
	abundance1 <- mglmOneGroup(y1+matrix(prior.count[j1],ntags,n1,byrow=TRUE),offset=offset.aug[j1],dispersion=dispersion)
	j2 <- group==pair[2]
	n2 <- sum(j2)
	if(n1==0) stop("No libraries for",pair[2])
	y2 <- y[,j2,drop=FALSE]
	abundance2 <- mglmOneGroup(y2+matrix(prior.count[j2],ntags,n2,byrow=TRUE),offset=offset.aug[j2],dispersion=dispersion)
	logFC <- (abundance2-abundance1)/log(2)

#	Equalize library sizes
	abundance <- mglmOneGroup(y,dispersion=dispersion,offset=offset)
	e <- exp(abundance)
	input.mean <- matrix(e,ntags,n1)
	output.mean <- input.mean*lib.size.average
	input.mean <- t(t(input.mean)*lib.size[j1])
	y1 <- q2qnbinom(y1,input.mean=input.mean,output.mean=output.mean,dispersion=dispersion)
	input.mean <- matrix(e,ntags,n2)
	output.mean <- input.mean*lib.size.average
	input.mean <- t(t(input.mean)*lib.size[j2])
	y2 <- q2qnbinom(y2,input.mean=input.mean,output.mean=output.mean,dispersion=dispersion)

	exact.pvals <- switch(rejection.region,
		doubletail=exactTestDoubleTail(y1,y2,dispersion=dispersion,big.count=big.count),
		deviance=exactTestByDeviance(y1,y2,dispersion=dispersion),
		smallp=exactTestBySmallP(y1,y2,dispersion=dispersion)
	)

	AveLogCPM <- object$AveLogCPM
	if(is.null(AveLogCPM)) AveLogCPM <- aveLogCPM(object)
	de.out <- data.frame(logFC=logFC, logCPM=AveLogCPM, PValue=exact.pvals)
	rn <- rownames(object$counts)
	if(!is.null(rn)) rownames(de.out) <- make.unique(rn)
	new("DGEExact",list(table=de.out, comparison=pair, genes=object$genes))
}
