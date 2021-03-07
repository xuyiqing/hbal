#' @title Double Selection
#' @aliases doubleSelection
#' @description Internal function called by \code{hbal} to perform double selection.
#' @param X                 covaraite matrix
#' @param W                 treatment indicator
#' @param Y                 outcome variable
#' @param grouping          groupings of covariates
#' @return resX, penalty.list, covar.keep
#' @author Yiqing Xu, Eddie Yang
#' @importFrom stats poly

doubleSelection <- function(
	X, 
	W, 
	Y,
	grouping
	){
	n <- 0
	penalty.list <- c()
	col.names <- colnames(X)
	resX <- getResidual(X, pos.list=grouping)
	start <- grouping[1] + 1
	end <- sum(grouping)
	P <- resX[,start:end]
	P.colnames <- col.names[start:end]
	# treatment
	cv.out <- cv.glmnet(P, W, alpha=1) 
	#lambda.min <- cv.out$lambda.min
	#t.coef <- coef(cv.out, s=lambda.min)
	t.coef <- coef(cv.out)
	t.coef <- t.coef@i 
	# outcome
	cv.out <- cv.glmnet(P, Y, alpha=1)
	#lambda.min <- cv.out$lambda.min
	#y.coef <- coef(cv.out, s=lambda.min)
	y.coef <- coef(cv.out)
	y.coef <- y.coef@i 
	all.coef <- sort(union(t.coef, y.coef))[-1]
	P <- P[,all.coef]
	P.colnames <- P.colnames[all.coef]
	P <- cbind(resX[,1:grouping[1]], P)
	colnames(P) <- c(col.names[1:grouping[1]], P.colnames)
	all.coef <- all.coef + grouping[1]
	for (i in 1:length(grouping)){
		nn <- n + grouping[i]
		penalty.list[i] <- sum(all.coef > n & all.coef <= nn)
		n <- nn
	}
	group.assignment <- penalty.list
	if(length(which(penalty.list==0))!=0){
		penalty.list <- penalty.list[-which(penalty.list==0)]
	}
	penalty.list <- c(grouping[1], penalty.list)
	covar.keep <- c(1:grouping[1], all.coef)
	out <- list(
		resX = P,
		penalty.list = penalty.list,
		covar.keep = covar.keep,
		select.group = group.assignment
		)
	return(out)
}