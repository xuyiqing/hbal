#' @title Estimating the ATT from an hbal object
#' @aliases att
#' @description \code{att} estimates the average treatment effect on the treated (ATT) from an 
#' hbal object returned by \code{hbal}. 
#' @usage att(hbalobject, method="lm_robust", ...)
#' @param hbalobject  an object of class \code{hbal} as returned by \code{hbal}.
#' @param method      estimation method for the ATT. Default is the Lin (2016) estimator. 
#' @param ...         arguments passed to lm_lin or  lm_robust
#' @details This is a wrapper for \code{lm_robust} and \code{lm_lin} from the \link{estimatr} package. 
#' @return A matrix of estimates with their robust standard errors
#' @importFrom estimatr lm_lin lm_robust
#' @importFrom stats as.formula
#' @author Yiqing Xu, Eddie Yang
#' @examples
#' #EXAMPLE 1
#' #Given W, Y, X, calculate the ATT through hbal
#' set.seed(92092)
#' W <- sample(0:1, 1000, replace=TRUE)
#' X <- data.frame(X1=rnorm(1000), X2=rnorm(1000))
#' Y <- rnorm(1000)
#' out <- hbal(Treatment = W, Y = Y, X = X)
#' lin <- att(out, method="lin")
#' @export

att <- function(
	hbalobject,
	method="lm_robust",
	...
	){
	if(class(hbalobject)!= "hbal" ){
     stop("hbalobject must be an hbal object from a call to hbal()")
    }
    w <- rep(NA, length(hbalobject$Treatment)) # get rid of "no visible binding for global variable" in check, see https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
	w[hbalobject$Treatment==1] <- 1
	w[hbalobject$Treatment==0] <- hbalobject$weights/sum(hbalobject$weights)*sum(hbalobject$Treatment==1)
	dat <- data.frame(hbalobject$mat, w=w, Y=hbalobject$Y, Treat=hbalobject$Treatment)
	if (method=="lin"){
		assignment <- as.formula(Y~Treat)
		covariates <- as.formula(paste0("~ ", paste0(colnames(hbalobject$mat), collapse=" + ")))
		out <- lm_lin(formula=assignment, covariates=covariates, weights=w, se_type="stata", data=dat, ...)
	}

	if (method=="lm_robust"){
		ff <- as.formula(paste0("Y~Treat+", paste0(colnames(hbalobject$mat), collapse=" + ")))
		out <- lm_robust(ff, weights=w, se_type="stata", data=dat, ...)
	}
	return(out)
}