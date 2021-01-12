#' @title Ridge Penalty Selection through Cross Validation
#' @aliases CrossValidate
#' @description Internal function called by \code{hbal} to select ridge penalties through cross-validation.
#' @param folds                   number of folds to perform cross validation. 
#' @param alpha                   alpha. Controls degree of regularization.
#' @param treatment               covariate matrix for treatment group.
#' @param fold.co                 fold assignments for control units.
#' @param fold.tr                 fold assignments for treated units.
#' @param coefs                   starting coefficients (lambda).
#' @param control                 covariate matrix for control group.
#' @param penalty                 list of hierarchical penalties.
#' @param constraint.tolerance    tolerance level for imbalance.
#' @param print.level             details of printed output.
#' @param base.weight             target weight distribution for the control units.
#' @param full.t                  (unresidualized) ovariate matrix for treatment group.
#' @param full.c                  (unresidualized) ovariate matrix for control group.
#' @param p                       positions to apply alpha
#' @param shuffle.treat           whether to create folds for the treated units
#' @return alpha, lambda
#' @importFrom stats na.omit coef
#' @importFrom glmnet cv.glmnet
#' @author Yiqing Xu, Eddie Yang

crossValidate <- function(
	folds=NULL,
	alpha=NULL,
	treatment = NULL,
	fold.co = NULL,
	fold.tr=NULL,
	coefs=NULL,
	control = NULL,
	penalty = NULL,
	constraint.tolerance = NULL,
	print.level = NULL,
	base.weight = NULL,
	full.t=NULL,
	full.c=NULL,
	p=NULL,
	shuffle.treat=NULL){

	res <- matrix(NA, length(alpha), folds) #store cross validation results
	rownames(res) <- as.character(1:length(alpha))
	coe <- NULL
	coef.store <- matrix(0, length(alpha), ncol(control))
	sub.coef <- matrix(c(coefs), folds, length(penalty), byrow=TRUE)

	# loop over each alpha value
	for (i in 1:length(alpha)){
	  penalty[p] <- alpha[i]
	  counter <- 0
		# loop over each fold for each alpha
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
			if(!is.null(sub.coef[k,]) && all(is.finite(sub.coef[k,]))) Coefs <- sub.coef[k,]

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

			if (class(out)!="try-error"){
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
					res[i,k] <- dif
					coef.store[i,] <- updateCoef(
						old.coef=coef.store[i,], 
						new.coef=coe, 
						counter=counter
						)
					counter <- counter + 1
				}
			}
			sub.coef[k,] <- coe
		}#end of inner loop
	}#end of outer loop
	# get rid of column (run) with all NA/NaN, unlucky partition
	col.keep <- !apply(!apply(res, 2, is.finite), 2, all)
	res <- res[,col.keep]
	res <-  na.omit(res)
	res.mean <- rowMeans(res)

	n <- as.numeric(rownames(res)[which(res.mean==min(res.mean))])[1]
	best.alpha <- alpha[n]
	lambda <- coef.store[n,]

	if (!all(is.finite(best.alpha))){
		best.alpha <- max(alpha)
		lambda <- coef.store[length(alpha),]
		warning("penalty replaced with max value")
	}

	return(
		list(
			best.alpha=best.alpha, 
			lambda=lambda
			)
		)
}



