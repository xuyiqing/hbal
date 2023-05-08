#' @title Ridge Penalty Selection through Cross Validation
#' @aliases CrossValidate
#' @description Internal function called by \code{hbal} to select ridge penalties through cross-validation.
#' @param group.alpha                   group.alpha. Controls degree of regularization.
#' @param penalty.pos             positions of user-supplied penalties.
#' @param penalty.val             values of user-supplied penalties.
#' @param group.exact             binary indicator of whether each covariate group should be penalized.
#' @param grouping                different groupings of the covariates.
#' @param folds                   number of folds to perform cross validation. 
#' @param treatment               covariate matrix for treatment group.
#' @param fold.co                 fold assignments for control units.
#' @param fold.tr                 fold assignments for treated units.
#' @param coefs                   starting coefficients (lambda).
#' @param control                 covariate matrix for control group.
#' @param constraint.tolerance    tolerance level for imbalance.
#' @param print.level             details of printed output.
#' @param base.weight             target weight distribution for the control units.
#' @param full.t                  (unresidualized) ovariate matrix for treatment group.
#' @param full.c                  (unresidualized) ovariate matrix for control group.
#' @param shuffle.treat           whether to create folds for the treated units
#' @return group.alpha, lambda
#' @importFrom stats na.omit coef
#' @importFrom glmnet cv.glmnet
#' @author Yiqing Xu, Eddie Yang

crossValidate <- function(
	group.alpha=NULL,
	penalty.pos=NULL,
	penalty.val=NULL,
	group.exact=NULL,
	grouping=NULL,
	folds=NULL,
	treatment = NULL,
	fold.co = NULL,
	fold.tr=NULL,
	coefs=NULL,
	control = NULL,
	constraint.tolerance = NULL,
	print.level = NULL,
	base.weight = NULL,
	full.t=NULL,
	full.c=NULL,
	shuffle.treat=NULL){

	if (any(!is.finite(group.alpha))){
		return(Inf)
	}
	
	res <- rep(NA, folds) #store cross validation results
	coe <- NULL

	# loop over each group.alpha value
	if (!is.null(group.exact)){
		group.alpha[which(group.exact==1)] <- 0
	}
	#penalty <- rep(c(0, group.alpha), times=grouping)
	penalty <- rep(group.alpha, times=grouping)
	if (!is.null(penalty.pos)){
		penalty[penalty.pos] <- penalty.val
	}

	sub.coef <- coefs
	Coefs <- coefs

	counter <- 0
	# loop over each fold for each group.alpha
	for (k in 1:folds){
		co.test.k <- which(fold.co==k)
		tr.test.k <- which(fold.tr==k)
		base.w <- base.weight[-co.test.k]
		train.control <- control[-co.test.k,]
		train.treat <- treatment[-tr.test.k,]
		train.total <- c(1, colMeans(train.treat))
		test.control <- control[co.test.k,]
		test.treat <- full.t[tr.test.k,]
		test.cc <- full.c[co.test.k,]
		if(!shuffle.treat){
			train.total <- c(1, colMeans(treatment))
			test.treat <- full.t
		}
		if(!is.null(sub.coef) && all(is.finite(sub.coef)) && max(abs(sub.coef))<=10) {
			Coefs <- sub.coef
		}	
		out <- try(
			hb(
				tr_total=as.matrix(train.total),
				co_x=train.control,
				coefs=as.matrix(Coefs),
				base_weight=as.matrix(base.w),
				alpha=as.matrix(penalty),   
				max_iterations=200,         
				constraint_tolerance=constraint.tolerance,
				print_level=print.level
				), 
			silent=TRUE
			)
		if (!inherits(out, "try-error")){
			coe <- out$coefs
			weights <- c(exp(test.control%*%coe))
			weights <- weights * base.weight[co.test.k]
			weights <- weights/sum(weights)
			test.weight <- all(is.finite(weights)) 
			test.coef <- all(is.finite(coe))
			if(test.weight && test.coef){
				test.cr.mean <- c(weights %*% test.cc)
				tr.t <- c(1, colMeans(test.treat))
				dif <- mean(abs(tr.t-test.cr.mean))
				res[k] <- dif
				counter <- counter + 1
			}
		}
		sub.coef <- coe
	}#end of inner loop

	oo <- ifelse(is.finite(mean(res)), mean(res), Inf)
	return(oo)
}



