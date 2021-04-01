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

	# treatment
	cv.out <- cv.glmnet(X, W, alpha=1) 
	t.coef <- coef(cv.out)
	t.coef <- t.coef@i 
	# outcome
	cv.out <- cv.glmnet(X, Y, alpha=1)
	y.coef <- coef(cv.out)
	y.coef <- y.coef@i 
	all.coef <- sort(union(t.coef, y.coef))[-1]
	all.coef <- sort(union(1:grouping[1], all.coef))

	X <- X[,all.coef]

	for (i in 1:length(grouping)){
		nn <- n + grouping[i]
		penalty.list[i] <- sum(all.coef > n & all.coef <= nn)
		n <- nn
	}

	if(length(which(penalty.list==0))!=0){
		penalty.list <- penalty.list[-which(penalty.list==0)]
	}

	out <- list(
		resX = X,
		penalty.list = penalty.list,
		covar.keep = all.coef
		)
	return(out)
}