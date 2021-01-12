#' @title Covariate Residualization
#' @aliases getResidual
#' @description Internal function called by \code{hbal} to residualize covariates.
#' @param P         covariate matrix
#' @param pos.list  groupings of covariates
#' @return Residualized covariates
#' @author Yiqing Xu, Eddie Yang

getResidual <- function(P, pos.list=NULL){
  n <- length(pos.list)-1
  for (i in 1:n){
    old <- sum(pos.list[1:i])
    new <- old + pos.list[i+1]
    residual.list <- apply(as.matrix(P[,(old+1):new]), 2, fLM, x=as.matrix(P[,1:old]))
    P[,(old+1):new] <- as.matrix(residual.list)
  }

  P <- as.matrix(P)
  return(P)
}