#' @title Hierarchically Regularized Entropy Balancing
#' @aliases hbal
#' @description \code{hbal} performs hierarchically regularized entropy balancing 
#' such that the covariate distributions of the control group match those of the 
#' treatment group. \code{hbal} automatically expands the covariate space to include
#' higher order terms and uses cross-validation to select variable penalties for the 
#' balancing conditions.
#' @usage hbal(Treatment, X, Y, base.weight = NULL, coefs = NULL ,
#'  max.iterations = 200, cv = TRUE, folds = 4, expand.degree = 3,
#'  ds = FALSE, alpha = NULL, constraint.tolerance = 1e-3, print.level = -1, 
#'  grouping = NULL, shuffle.treat=TRUE, exclude=NULL)
#' @param Treatment            a numeric, binary vector of treatment status. 1 should denote treatment observations and 0 should denote control observations.
#' @param X                    matrix of covariates. When expand = \code{TRUE}, X is serially expanded to include higher-order terms of X.
#' @param Y                    numeric vecor of outcome variable. Used only by subsequent estimation. See \code{att}.
#' @param base.weight          target weight distribution for the control units.
#' @param coefs                initial coefficients for the reweighting algorithm (lambdas).
#' @param max.iterations       maximum number of iterations. Default is 200.
#' @param cv                   whether to use cross validation. Default is \code{TRUE}.
#' @param folds                number of folds for cross validation. Only used when cv is \code{TRUE}.
#' @param expand.degree        degree of series expansion. 0 means no expansion. Default is 3.
#' @param ds                   whether to perform double selection prior to balancing. Default is \code{FALSE}.
#' @param alpha                vector of ridge penalties used in grid search during cross-validation. Only used when cv is \code{TRUE}.
#' @param constraint.tolerance tolerance level for overall imbalance. Default is 1e-3.
#' @param print.level          details of printed output.
#' @param grouping             different groupings of the covariates. Must be specified if expand is \code{FALSE}.
#' @param shuffle.treat        whether to use cross-validation on the treated units. Default is \code{TRUE}.
#' @param exclude              list of covariate name pairs or triplets to be excluded.
#' @details In the simplest set-up, user can just pass in \{Treatment, X, Y\}. The default settings will serially expand
#' X to include higher order terms, hierarchically residualize these terms, perform double selection to only keep the relevant
#' variables and use cross-validation to select penalities for different groupings of the covariates. 
#' @return 
#' An list object of class \code{hbal} with the following elements:
#' \item{coefs}{vector that contains coefficients from the reweighting algorithm.}
#' \item{mat}{matrix of serially expanded covariates if expand=\code{TRUE}. Otherwise, the original covariate matrix is returned.}
#' \item{penalty}{vector of ridge penalties used for each covariate} 
#' \item{weights}{vector that contains the control group weights assigned by hbal.}
#' \item{W}{vector of treatment status}
#' \item{Y}{vector of outcome}
#' @author Yiqing Xu, Eddie Yang
#' @importFrom stats var
#' @useDynLib hbal, .registration = TRUE
#' @examples
#' # Example 1
#' set.seed(92092)
#' N <- 500
#' X1 <- rnorm(N)
#' X2 <- rbinom(N,size=1,prob=.5)
#' X <- cbind(X1, X2)
#' treat <- rbinom(N, 1, prob=0.5) # Treatment indicator
#' y <- X[,1] + X[,2] + rnorm(N) # Outcome
#' out <- hbal(Treatment = treat, Y = y, X = X)
#' summary(hbal::att(out))
#' 
#' # Example 2
#' ## Simulation from Kang and Shafer (2007).
#' library(MASS)
#' set.seed(92092)
#' n <- 500
#' X <- mvrnorm(n, mu = rep(0, 4), Sigma = diag(4))
#' prop <- 1 / (1 + exp(X[,1] - 0.5 * X[,2] + 0.25*X[,3] + 0.1 * X[,4]))
#' # Treatment indicator
#' treat <- rbinom(n, 1, prop)
#' # Outcome
#' y <- 210 + 27.4*X[,1] + 13.7*X[,2] + 13.7*X[,3] + 13.7*X[,4] + rnorm(n)
#' # Observed covariates
#' X.mis <- cbind(exp(X[,1]/2), X[,2]*(1+exp(X[,1]))^(-1)+10, 
#'     (X[,1]*X[,3]/25+.6)^3, (X[,2]+X[,4]+20)^2)
#' out <- hbal(Treatment = treat, Y = y, X = X.mis)
#' summary(att(out))
#' @export

hbal <- function(
	Treatment,
	X,
	Y,
	base.weight = NULL,
	coefs = NULL ,
	max.iterations = 200,
	cv=TRUE,
	folds=4,
	expand.degree=3,
	ds=FALSE,
	alpha=NULL,
	constraint.tolerance = 1e-3,
	print.level=-1,
	grouping=NULL,
	shuffle.treat=TRUE,
	exclude=NULL
	){

	# ntreated: number of treated units
	# ncontrols: number of control units
	# full: standardized design matrix (X) to be used for cross validation
	# full.t: standardized covariates for the treated units
	# full.c: standardized covariates for the control units, also has a column of 1s for normalizing constraint
	# co.x: covariates for the control group, also has a column of 1s
	# tr.total: mean of covariates of the treatment group, first element (1) for the normalizing constraint
	# hyperpara: \alpha_{k}, multiplier on the regularization term
	# fold.num.co: sequence of 1:folds of length ceiling(ncontrols/folds) to be sampled from, control group
	# fold.co: assign a fold number to each control unit
	# fold.num.tr: sequence of 1:folds of length ceiling(ntreated/folds) to be sampled from, treatment group
	# fold.tr: assign a fold number to each treated unit
	# lambda: lambda in the paper, Lagrangian multipliers, will be used to store the returned Lagrangian multipliers 
	# penalty: sequence of \alpha of length equal to the length of coefs
	# n.group: number of higher order groups 
	# end, start, p: end and starting points to search for \alpha
	# weights: solution weights
	# cc: used to store final Lagrangian multipliers 

	# set up elements
	mcall <- match.call()
	X  <- as.matrix(X)
	if (is.null(colnames(X))) colnames(X) <- paste0("X", 1:ncol(X))
	if (sum(is.na(X))!=0) stop("data contain missing values")
	if (is.null(alpha) & cv==TRUE) alpha <- c(0,exp(seq(log(0.01), log(1000), length.out = 24))) # 25 alpha values distributed exponentially (0, 100) for grid search, this controls the degree of regularization for each group
	if (!is.numeric(Treatment)) stop("Treatment indicator needs to be numeric")
	ntreated  <- sum(Treatment==1)
	ncontrols <- sum(Treatment==0)  
	
	if (expand.degree > 0) {# series expansion of the covariates
		# need to residualize here actually so that (expand=TRUE, ds=FALSE) also residualize the covariates
		expand <- covarExpand(X, exp.degree=expand.degree, treatment=Treatment, exclude=exclude)
		X <- scale(expand$mat)
		full <- X # need to keep a copy of the standardized X for cross-validation. Normalization so that means and variances of covariates on the same scale 
		full.t <- full[Treatment==1,]
		full.c <- full[Treatment==0,]
		grouping <- expand$grouping
		#if (!ds) X <- scale(getResidual(X, pos.list=grouping))
	} else{
		X <- scale(X)
		full <- X
		full.t <- full[Treatment==1,]
		full.c <- full[Treatment==0,]
	}

	if (ds){
		selected <- doubleSelection(X=X, W=Treatment, Y=Y, grouping=grouping)
		X <- scale(selected$resX)
		grouping <- selected$penalty.list
		full <- full[,selected$covar.keep]
		full.t <- full.t[,selected$covar.keep]
		full.c <- full.c[,selected$covar.keep]
	}
	full.c <- cbind(1, full.c)

	if(is.null(grouping)) grouping <- ncol(X)
	co.x <- X[Treatment==0,] # control group
	co.x <- cbind(rep(1,ncontrols),co.x)
	tr.total <- c(1, colMeans(X[Treatment==1,,drop=FALSE]))
	grouping[1] <- grouping[1] + 1 # need one more for normalizing constraint

	# Checks
	if (sum(grouping) != ncol(co.x)) stop("Sum of grouping must be equal to ncol(X)")
	if (sum(Treatment != 1 & Treatment != 0) > 0) stop("Treatment indicator ('Treatment') must be a logical variable, TRUE (1) or FALSE (0)")
	if (var(Treatment) == 0) stop("Treatment indicator ('Treatment') must contain both treatment and control observations")
	if (sum(is.na(Treatment)) > 0) stop("Treatment contains missing data")
	if (length(Treatment) != nrow(X)) stop("length(Treatment) != nrow(X)")
	if (sum(is.na(X)) > 0) stop("X contains missing data")
	if (length(max.iterations) != 1 ) stop("length(max.iterations) != 1")
	if (length(constraint.tolerance) != 1 ) stop("length(constraint.tolerance) != 1")
	if (is.null(base.weight)) base.weight = rep(1, ncontrols)
	if (length(base.weight) !=  ncontrols) stop("length(base.weight) !=  number of controls  sum(Treatment==0)")
	if (qr(co.x)$rank != ncol(co.x)) stop("collinearity in covariate matrix for controls (remove collinear covariates)")
	if (is.null(coefs)) coefs = c(log(tr.total[1]/sum(base.weight)),rep(0,(ncol(co.x)-1)))
	if (length(coefs) != ncol(co.x)) stop("coefs needs to have same length as number of covariates plus one")
	if (is.null(alpha) & cv==FALSE) alpha <- c(0,length(coefs))
	if (is.logical(cv)==FALSE) stop("cv needs to be of type logical")
	if (cv==TRUE && length(grouping)==1){
		cv <- FALSE
		warning("length(grouping)==1, reverting to cv = FALSE. \nEither double selection selected 0 higher order term or the supplied grouping has length 1")
	}
	if (print.level >= 1){
	  cat("Data Setup\nCovariate Adjustment:", colnames(X), "\n" )
	  cat("\n")
	  cat("Control to Treatment ratio = ", ncontrols/ntreated, "\n")
	}

	# Main
	if (cv==FALSE){
		z <- hb(
			tr_total=as.matrix(tr.total),
			co_x=co.x,
			coefs=as.matrix(coefs),
			base_weight=as.matrix(base.weight),
			alpha=alpha, 
			max_iterations=max.iterations,
			constraint_tolerance=constraint.tolerance,
			print_level=print.level
			)
		hyperpara <- list(hyperpara = NULL)
	} else{
		fold.num.co <- rep(1:folds, ceiling(ncontrols/folds))
		fold.co <- sample(fold.num.co, ncontrols, replace=F) 
		fold.num.tr <- rep(1:folds, ceiling(ntreated/folds))
		fold.tr <- sample(fold.num.tr, ntreated, replace=F)

		#initilize penalty vector to max(alpha), except for penalty for linear terms
		lambda <- coefs
		penalty <- rep(max(alpha), length(coefs))
		penalty[1:grouping[1]] <- 0
		n.group <- length(grouping) - 1 # no. of higher order groups

		for (i in 1:n.group){
			i <- i + 1 # not regularizing linear terms
			end <- sum(grouping[1:i])
			start <- end - grouping[i] +1
			p <- start:end #positions to apply penalty
			res <- crossValidate(
				folds=folds,
				alpha=alpha,
				fold.co=fold.co,
				fold.tr=fold.tr,
				treatment = X[Treatment==1,],
				control=co.x,
				coefs = lambda,
				penalty=penalty,
				base.weight=base.weight,
				constraint.tolerance=constraint.tolerance,
				print.level=print.level,
				full.t=full.t,
				full.c=full.c,
				p=p,
				shuffle.treat=shuffle.treat
				)

			lambda <- res$lambda
			penalty[p] <- res$best.alpha
		} #end of for loop
	} #end of if, else

	if (cv==FALSE){
		weights<-z$Weights_ebal
		cc <- z$coefs
	}else{
		weights <- c(exp(co.x%*%res$lambda))
		weights <- weights * base.weight
		weights <- weights/sum(weights)
		cc <- res$lambda
	}
	if (expand.degree>1 & ds){
		out.mat <- expand$mat[,selected$covar.keep]
	}
	if (expand.degree>1 & !ds){
		out.mat <- expand$mat
	}
	if (expand.degree<1 & !ds){
		out.mat <- X
	}
	if (ds){
		group.assignment <- selected$select.group
	}else{
		grouping[1] <- grouping[1]-1
		group.assignment <- grouping
	}
	out <- list(weights=weights, 
				coefs=cc, 
				Treatment=Treatment, 
				Y=Y, 
				mat=out.mat, 
				group.assignment=group.assignment,
				penalty=NULL,
				call=mcall)

	if (cv==TRUE){
		out[["penalty"]] <- penalty[-1]
	}

	class(out) <- "hbal"

	return(out)
}