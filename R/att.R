#' @title Estimating the ATT from an hbal object
#' @aliases att
#' @description \code{att} estimates the average treatment effect on the treated (ATT) from an 
#' hbal object returned by \code{hbal}. 
#' @usage att(hbalobject, method="lm_robust", dr=TRUE, ...)
#' @param hbalobject  an object of class \code{hbal} as returned by \code{hbal}.
#' @param method      estimation method for the ATT. Default is the Lin (2016) estimator. 
#' @param dr      	  doubly robust, whether an outcome model is included in estimating the ATT.
#' @param ...         arguments passed to lm_lin or lm_robust
#' @details This is a wrapper for \code{lm_robust} and \code{lm_lin} from the \link{estimatr} package. 
#' @return A matrix of estimates with their robust standard errors
#' @importFrom estimatr lm_lin lm_robust
#' @importFrom stats as.formula
#' @author Yiqing Xu, Eddie Yang
#' @examples
#' #EXAMPLE 1
#' set.seed(1984)
#' N <- 500
#' X1 <- rnorm(N)
#' X2 <- rbinom(N,size=1,prob=.5)
#' X <- cbind(X1, X2)
#' treat <- rbinom(N, 1, prob=0.5) # Treatment indicator
#' y <- 0.5 * treat + X[,1] + X[,2] + rnorm(N) # Outcome
#' dat <- data.frame(treat=treat, X, Y=y)
#' out <- hbal(Treat = 'treat', X = c('X1', 'X2'), Y = 'Y', data=dat)
#' sout <- summary(att(out))
#' @export

att <- function(
	hbalobject,
	method="lm_robust",
	dr=TRUE,
	...
	){
	if(class(hbalobject)!= "hbal" ){
     stop("hbalobject must be an hbal object from a call to hbal()")
    }
    elpss <- list(...)
    
    w <- rep(NA, length(hbalobject$Treatment)) # get rid of "no visible binding for global variable" in check, see https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
	w[hbalobject$Treatment==1] <- 1
	w[hbalobject$Treatment==0] <- hbalobject$weights/sum(hbalobject$weights)*sum(hbalobject$Treatment==1)
	dat <- data.frame(hbalobject$mat, w=w, Y=hbalobject$Y, Treat=hbalobject$Treatment)
	if (method=="lin"){
		covariates <- as.formula(paste0("~ ", paste0(colnames(hbalobject$mat), collapse=" + ")))
		if (is.null(elpss[['se_type']])){
			out <- lm_lin(formula=as.formula(Y~Treat), covariates=covariates, weights=w, data=dat, se_type='stata', ...)
		}else{
			out <- lm_lin(formula=as.formula(Y~Treat), covariates=covariates, weights=w, data=dat, ...)
		}
	}

	if (method=="lm_robust"){
		if(dr){
			ff <- as.formula(paste0("Y~Treat+", paste0(colnames(hbalobject$mat), collapse=" + ")))
		}else{
			ff <- as.formula("Y~Treat")
		}

		if (is.null(elpss[['se_type']])){
			out <- lm_robust(ff, data=dat, weights=w, se_type='stata', ...)
		}else{
			out <- lm_robust(ff, data=dat, weights=w, ...)
		}
	}
	return(out)
}