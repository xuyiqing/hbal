#' @title Residualization
#' @importFrom RcppEigen fastLm
#' @author Yiqing Xu, Eddie Yang
#' @examples
#' @export

fLM <- function(y, x){
	out <- RcppEigen::fastLm(X=x, y=y)
	return(out$residuals)
}


getResidual <- function(P, pos.list=NULL){
	n <- length(pos.list)-1
  for (i in 1:n){
    old <- sum(pos.list[1:i])
    new <- old + pos.list[i+1]
    residual_list <- apply(as.matrix(P[,(old+1):new]), 2, fLM, x=as.matrix(P[,1:old]))
    P[,(old+1):new] <- as.matrix(residual_list)
  }
	#P <- scale(P)
	P <- as.data.frame(P)
	return(P)
}








