#' @title Estimating the ATT from an hbal object
#' @aliases att
#' @description \code{att} estimates the average treatment effect on the treated (ATT) from an 
#' hbal object returned by \code{hbal}. 
#' @usage att(hbalobject, method="lm_robust", dr=TRUE, displayAll=FALSE, ...)
#' @param hbalobject  an object of class \code{hbal} as returned by \code{hbal}.
#' @param method      estimation method for the ATT. Default is the Lin (2016) estimator. 
#' @param dr      	  doubly robust, whether an outcome model is included in estimating the ATT.
#' @param displayAll  only displays treatment effect by default.
#' @param ...         arguments passed to lm_lin or lm_robust
#' @details This is a wrapper for \code{lm_robust} and \code{lm_lin} from the \link{estimatr} package. 
#' @return A matrix of estimates with their robust standard errors
#' @importFrom estimatr lm_lin lm_robust
#' @importFrom stats as.formula
#' @importFrom generics tidy
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
	displayAll=FALSE,
	...
	){
	if(!inherits(hbalobject, "hbal")){
		stop("hbalobject must be an hbal object from a call to hbal()")
    }
    elpss <- list(...)
	if (dr == FALSE & method == "lin") {
		method <- "lm_robust"
	}
	if (is.null(hbalobject$Y)==TRUE) {
		stop("The outcome variable is missing in hbalobject; please add Y when calling hbal()")
	}
	# extract data 
	mat <-   hbalobject$mat 
    Y <- hbalobject$Y # name of the outcome variable
	Tr <- hbalobject$Treat # name of the treatment variable
	Covar <- colnames(mat)
	dat <- cbind.data.frame(hbalobject$Outcome, hbalobject$Treatment, mat)
	colnames(dat) <- c(Y, Tr, Covar)
	w <- "w" # this line is useless; just to get around CRAN checker
	dat$w <- hbalobject$weights
		# linear regression
	if (method=="lm_robust"){
		if(dr){
			ff <- as.formula(paste0(Y, ' ~ ', Tr, ' + ', paste0(Covar, collapse=" + ")))
		}else{
			ff <- as.formula(paste0(Y, ' ~ ', Tr))
		}

		if (is.null(elpss[['se_type']])){
			out <- lm_robust(ff, data=dat, weights=w, se_type='stata', ...)
		}else{
			out <- lm_robust(ff, data=dat, weights=w, ...)
		}
	}
	# Lin (2013): regression with interactions
	if (method=="lm_lin"){
		covariates <- as.formula(paste0("~ ", paste0(Covar, collapse=" + ")))
		if (is.null(elpss[['se_type']])){
			out <- lm_lin(formula = as.formula(paste0(Y, ' ~ ', Tr)), covariates=covariates, weights=w, data=dat, se_type='stata', ...)
		}else{
			out <- lm_lin(formula = as.formula(paste0(Y, ' ~ ', Tr)), covariates=covariates, weights=w, data=dat, ...)
		}
	}
	if (displayAll == FALSE)
	{
		out <- tidy(out)[2, -c(1, 9), drop = FALSE]
		colnames(out) <-  c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "CI Lower", "CI Upper", "DF")
		rownames(out) <- Tr	
		return(out)	  
	}
	else
	{
	  return(out)
	}
	
}