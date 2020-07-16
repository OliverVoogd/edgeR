estimateGLMRobustDisp <-
function (y, design = NULL, prior.df = 10, update.trend = TRUE, trend.method = "bin.loess", maxit = 6, k = 1.345, residual.type = "pearson", verbose = FALSE, record = FALSE)
{

	if (!is(y, "DGEList")) 
	stop("Input must be a DGEList.")
	y$weights <- array(1, dim(y))
	if(is.null(y$trended.dispersion)) y <- estimateGLMTrendedDisp(y,design = design, method = trend.method)
	if(is.null(y$tagwise.dispersion)) y <- estimateGLMTagwiseDisp(y,design = design, prior.df = prior.df)   
	if(record) y <- .getRecord(y, i = 0, weights = y$weights)
	for (i in seq_len(maxit)){
		if (verbose) 
		message(paste0("Iteration ", i, ": Re-fitting GLM. "), 
				appendLF = FALSE)
		fit <- glmFit(y, design = design)
		res <- .calcResid(fit, residual.type = residual.type)
		y$weights <- .psi.huber.matrix(res, k = k)
		y$AveLogCPM <- aveLogCPM(y,dispersion=y$trended.dispersion)
		if(update.trend){
			if (verbose) 
			message("Re-estimating trended dispersion.")	
			y <- estimateGLMTrendedDisp(y, design = design, method = trend.method)
		}
		if (verbose) 
		message("Re-estimating tagwise dispersion.")
		y <- estimateGLMTagwiseDisp(y, design = design, prior.df=prior.df)
		if(record) y <- .getRecord(y, i = i, res = res, weights = y$weights, fit = fit)
		
	}	
	y	
	
}




.calcResid <- function(f,residual.type=c("pearson", "anscombe","deviance")) 
{	residual.type <- match.arg(residual.type)
	resAns <- function(y,mu,disper)
	{ 
		res <- function(y,mu,disper)    
		{     
			f <- function(x,disper) {(x*(1+disper*x))^(-1/3)}
			const <- (f(mu,disper=disper))^(1/2)
			if(mu==0)
			  out <- 0 
			else 
			  out <- const*integrate(f,mu,y,disper=disper)$value 
			out
		}
		resV<- Vectorize(res,vectorize.args=c("y","mu","disper"))
		out <- matrix(resV(y,mu,disper),nrow=nrow(y))
	}	
	
	resDev <- function(y,mu,disper)
	{ 
		y <- y+1e-5
		r <- 2*(y*log(y/mu)+(y+1/disper)*log((mu+1/disper)/(y+1/disper)))
        #r[y==0] <- 0
		r[mu==0] <- 0
		sign(y-mu)*sqrt(r)
	}
	
	mu <- f$fitted.values
	disp <- expandAsMatrix(f$dispersion,dim(mu))
	yi <- f$counts
	res <- switch(residual.type,anscombe=resAns(yi,mu,disp),
					  pearson = {(yi - mu)/sqrt((mu * (1+(disp)*mu)))},
					  deviance = resDev(yi,mu,disp))
	res[mu==0] <- 0
	res
}



.psi.huber.matrix <- function (u,k=1.345)
{
	z <- k/abs(u)
	z[abs(u) <= k] <- 1
	z
}



.getRecord <-
function(y, i, res = NULL, weights = NULL, fit = NULL)
{
	iteration <- paste0("iteration_", i)
	if(is.null(y$record))	record <- list()
	else record <- y$record	
	if(!is.null(y$AveLogCPM)) record$AveLogCPM[[iteration]] <- y$AveLogCPM			 
	if(!is.null(y$trended.dispersion)) record$trended.dispersion[[iteration]] <- y$trended.dispersion
	if(!is.null(y$tagwise.dispersion)) record$tagwise.dispersion[[iteration]] <- y$tagwise.dispersion
	if(!is.null(weights)) record$weights[[iteration]] <- weights
	if(!is.null(res)) record$res[[iteration]] <- res
	if(!is.null(fit)) record$mu[[iteration]] <- fit$fitted.values
	y$record <- record
	y		
}