#' @title Update lambda
#' @aliases updateCoef
#' @description Internal function called by \code{hbal} to residualize covariates.
#' @param old.coef previous coefficients
#' @param new.coef new coefficients
#' @param counter  which fold in CV
#' @return updated coefficients
#' @author Yiqing Xu, Eddie Yang

updateCoef <- function(old.coef, new.coef, counter){
	updated.coef <- (old.coef * counter + new.coef)/(counter+1)
	return(updated.coef)
}