#' @title Covariate Residualization
#' @aliases fLM
#' @description Wrapper for RcppEigen's fastLm. Internal function called by \code{hbal} to residualize covariates.
#' @param y x in RcppEigen::fastLm
#' @param x y in RcppEigen::fastLm
#' @return Residualized covariates
#' @author Yiqing Xu, Eddie Yang
#' @importFrom RcppEigen fastLmPure

fLM <- function(y, x){
	out <- RcppEigen::fastLmPure(X=x, y=y, method = 1L)
	return(out$residuals)
}