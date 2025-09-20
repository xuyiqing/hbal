#' @title Estimating the ATT from an hbal object
#' @aliases att
#' @description \code{att} estimates the average treatment effect on the treated (ATT) from an 
#' hbal object returned by \code{hbal}. 
#' @usage att(hbalobject, method="lm_robust", dr=TRUE, displayAll=FALSE, alpha=0.9, ...)
#' @param hbalobject  an object of class \code{hbal} as returned by \code{hbal}.
#' @param method      estimation method for the ATT. Default is the Lin (2016) estimator. 
#' @param dr      	  doubly robust, whether an outcome model is included in estimating the ATT.
#' @param displayAll  only displays treatment effect by default.
#' @param alpha       tuning paramter for glmnet
#' @param ...         arguments passed to lm_lin or lm_robust
#' @details This is a wrapper for \code{lm_robust} and \code{lm_lin} from the \link{estimatr} package. 
#' @return A matrix of estimates with their robust standard errors
#' @importFrom estimatr lm_lin lm_robust
#' @importFrom stats as.formula
#' @importFrom stats predict
#' @importFrom generics tidy
#' @importFrom glmnet cv.glmnet
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
	alpha = 0.9,
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
	# Athey (2018): Approximate residual balancing
	if (method == "elnet"){
	  # Estimate est0
	  XW <- as.matrix(dat[hbalobject$Treatment==0, Covar, drop = FALSE])
	  XT <- as.matrix(dat[hbalobject$Treatment==1, Covar, drop = FALSE])
	  YW <- as.matrix(dat[hbalobject$Treatment==0, Y, drop = FALSE])
	  YT <- as.matrix(dat[hbalobject$Treatment==1, Y, drop = FALSE])
	  lasso.fit <-  glmnet::cv.glmnet(XW, YW, alpha = alpha)
	  mu.lasso <-  predict(lasso.fit, newx = matrix(colMeans(XT), 1, ncol(XT)))
	  residuals <-  YW - predict(lasso.fit, newx = XW)
	  gamma <- hbalobject$weights.co/sum(hbalobject$weights.co)
	  mu.residual <- sum(gamma * residuals)
	  ## degrees of freedom correction
	  var.hat <-  sum(gamma^2 * residuals^2) * length(gamma) / max(1, length(gamma) - sum(coef(lasso.fit) != 0))
	  mu.hat <-  mu.lasso + mu.residual
	  est0 <- c(mu.hat, var.hat)
	  df0 <- max(1, length(gamma) - sum(coef(lasso.fit) != 0))
	  # Estimate est1
	  lasso.fit <-  glmnet::cv.glmnet(XT, YT, alpha = alpha)
	  residuals <-  YT - predict(lasso.fit, newx = XT)
	  var.hat <- mean(residuals^2) / max(1, length(YT) - sum(coef(lasso.fit) != 0))
	  est1 <-  c(mean(YT), var.hat)
	  df1 <- max(1, length(YT) - sum(coef(lasso.fit) != 0))
	  # tau and var
	  tau.hat <-  est1[1] - est0[1]
	  var.hat <-  est1[2] + est0[2]
	  df <- df0 + df1
	  # out
	  out <- data.frame("Estimate" = tau.hat, "Std. Error" = sqrt(var.hat), "t value" = tau.hat/sqrt(var.hat), "Pr(>|t|)" = 2*stats::pt(-tau.hat/sqrt(var.hat), df = df), "CI Lower" = tau.hat - 1.96 * sqrt(var.hat), "CI Upper" = tau.hat + 1.96 * sqrt(var.hat), "DF" = df)
	  colnames(out) <-  c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "CI Lower", "CI Upper", "DF")
	  rownames(out) <- Tr	
	  return(out)
	}
	# DisplayAll Option of ATT
	if (displayAll == FALSE){
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