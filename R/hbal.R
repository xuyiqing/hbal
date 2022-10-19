#' @title Hierarchically Regularized Entropy Balancing
#' @aliases hbal
#' @description \code{hbal} performs hierarchically regularized entropy balancing 
#' such that the covariate distributions of the control group match those of the 
#' treatment group. \code{hbal} automatically expands the covariate space to include
#' higher order terms and uses cross-validation to select variable penalties for the 
#' balancing conditions.
#' @usage hbal(data, Treat, X, Y, base.weight = NULL, coefs = NULL ,
#'      max.iterations = 200, cv = TRUE, folds = 4, 
#'      expand.degree = 3, ds = FALSE, alpha = NULL, 
#'      group.alpha=NULL, constraint.tolerance = 1e-3, 
#'      print.level = -1, grouping = NULL, 
#'      shuffle.treat=TRUE, exclude=NULL, seed=NULL)
#' @param data                 a dataframe that contains the treatment, outcome, and covariates.   
#' @param Treat                a character string of the treatment variable.
#' @param X                    a character vector of covariate names to balance on.
#' @param Y                    a character string of the outcome variable.
#' @param base.weight          target weight distribution for the control units.
#' @param coefs                initial coefficients for the reweighting algorithm (lambdas).
#' @param max.iterations       maximum number of iterations. Default is 200.
#' @param cv                   whether to use cross validation. Default is \code{TRUE}.
#' @param folds                number of folds for cross validation. Only used when cv is \code{TRUE}.
#' @param expand.degree        degree of series expansion. 0 means no expansion. Default is 3.
#' @param ds                   whether to perform double selection prior to balancing. Default is \code{FALSE}.
#' @param alpha                named vector of ridge penalties.
#' @param group.alpha          binary indicator of whether each covariate group should be penalized.
#' @param constraint.tolerance tolerance level for overall imbalance. Default is 1e-3.
#' @param print.level          details of printed output.
#' @param grouping             different groupings of the covariates. Must be specified if expand is \code{FALSE}.
#' @param shuffle.treat        whether to use cross-validation on the treated units. Default is \code{TRUE}.
#' @param exclude              list of covariate name pairs or triplets to be excluded.
#' @param seed                 random seed to be set. Set random seed when cv=\code{TRUE} for reproducibility.
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
#' @importFrom nloptr nloptr
#' @useDynLib hbal, .registration = TRUE
#' @references Xu, Y., & Yang, E. (2022). Hierarchically Regularized Entropy Balancing. Political Analysis, 1-8. doi:10.1017/pan.2022.12
#' @examples
#' # Example 1
#' set.seed(1984)
#' N <- 500
#' X1 <- rnorm(N)
#' X2 <- rbinom(N,size=1,prob=.5)
#' X <- cbind(X1, X2)
#' treat <- rbinom(N, 1, prob=0.5) # Treatment indicator
#' y <- 0.5 * treat + X[,1] + X[,2] + rnorm(N) # Outcome
#' dat <- data.frame(treat=treat, X, Y=y)
#' out <- hbal(Treat = 'treat', X = c('X1', 'X2'), Y = 'Y', data=dat)
#' summary(hbal::att(out))
#' 
#' # Example 2
#' ## Simulation from Kang and Shafer (2007).
#' library(MASS)
#' set.seed(1984)
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
#' dat <- data.frame(treat=treat, X.mis, Y=y)
#' out <- hbal(Treat = 'treat', X = c('X1', 'X2', 'X3', 'X4'), Y='Y', data=dat)
#' summary(att(out))
#' @export

hbal <- function(
	data,
	Treat,
	X,
	Y,
	base.weight=NULL,
	coefs=NULL ,
	max.iterations=200,
	cv=TRUE,
	folds=4,
	expand.degree=3,
	ds=FALSE,
	alpha=NULL,
	group.alpha=NULL,
	constraint.tolerance=1e-3,
	print.level=-1,
	grouping=NULL,
	shuffle.treat=TRUE,
	exclude=NULL,
	seed=NULL
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
	if(!is.null(seed)) set.seed(seed)

	var.names <- colnames(data)
	if(!all(c(Treat, X, Y) %in% var.names)) stop("Some variable(s) specified are not in the data")

	num_col <- colSums(sapply(data[,X], is.finite))==nrow(dat)
	if (sum(num_col)!=length(num_col)) warning(paste0('Variables: ', X[!num_col], ' are dropped because they are not numeric/finite'))
	X <- X[num_col]

	X  <- raw <- as.matrix(data[,X])

	if (is.null(colnames(X))) colnames(X) <- paste0("X", 1:ncol(X))
	if (sum(is.na(X))!=0) stop("data contain missing values")

	Treatment <- unlist(data[,Treat])
	if (!is.numeric(Treatment)) stop("Treatment variable needs to be numeric")
	if (var(Treatment)==0) stop("Treatment variable only contains treated/control units")
	ntreated  <- sum(Treatment==1)
	ncontrols <- sum(Treatment==0)  
	
	if(is.null(grouping)) grouping <- ncol(X)
	if(length(expand.degree)!=1) stop("expand.degree needs to be of length 1")
	if (!expand.degree %in% c(1, 2, 3)) stop("expand.degree needs to be one of c(1, 2, 3)")
	if (expand.degree > 1) {# series expansion of the covariates
		expand <- covarExpand(X, exp.degree=expand.degree, treatment=Treatment, exclude=exclude)
		grouping <- expand$grouping
		if (0 %in% grouping){
    		grouping <- grouping[-which(grouping==0)]
		}
		X <- scale(expand$mat)
	} else{
		X <- scale(X)
	}

	full <- X # need to keep a copy of the standardized X for cross-validation. Normalization so that means and variances of covariates on the same scale 
	full.t <- full[Treatment==1,]
	full.c <- full[Treatment==0,]

	if (ds){
		if (expand.degree > 1) {X.ds <- expand$mat} else {X.ds <- raw}
		selected <- doubleSelection(X=X.ds, W=Treatment, Y=unlist(data[,Y]), grouping=grouping)
		X <- scale(selected$resX)
		grouping <- selected$penalty.list
	}
	
	if (length(grouping)==1){
		cv <- FALSE
		if (print.level>0){
			cat("length(grouping)==1, setting cv=FALSE")
		}
	}

	full.c <- cbind(rep(1,ncontrols), full.c)

	co.x <- X[Treatment==0,] # control group
	co.x <- cbind(rep(1,ncontrols),co.x)
	tr.total <- c(1, colMeans(X[Treatment==1,,drop=FALSE]))
	grouping[1] <- grouping[1] + 1 # need one more for normalizing constraint

	if(!is.null(group.alpha)){
		group.alpha <- group.alpha[-1]
	}
	if (!is.null(alpha)){
		penalty.names <- names(alpha)
		penalty.val <- unname(alpha)
		penalty.pos <- match(penalty.names, colnames(X))+1
	}else{
		penalty.names <- penalty.val <- penalty.pos <- NULL
	}
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
	if (is.null(alpha) & cv==FALSE) alpha <- rep(0,length(coefs))
	if (is.logical(cv)==FALSE) stop("cv needs to be of type logical")
	if (cv==TRUE && length(grouping)==1){
		cv <- FALSE
		warning("length(grouping)==1, reverting to cv = FALSE. \nEither double selection selected 0 higher order term or the supplied grouping has length 1")
	}
	if (print.level >= 1){
		cat("Data Setup\nCovariate Adjustment:", colnames(X), "\n\n" )
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


		min.c <- nloptr(x0 = rep(1, length(grouping)-1),
                eval_f = crossValidate,
                lb = rep(0, length(grouping)-1),
                ub = rep(100, length(grouping)-1),
                opts = list('algorithm'='NLOPT_LN_COBYLA',
                            'maxeval' =200,
			    'ftol_rel'=1e-3,
			    'ftol_abs'=1e-5,
			    'xtol_abs'=1e-3,
                            'print_level'=print.level),
                penalty.pos=penalty.pos,
                penalty.val=penalty.val,
                group.alpha=group.alpha,
                grouping=grouping,
                folds=folds,
                treatment = X[Treatment==1,],
                fold.co = fold.co,
                fold.tr=fold.tr,
                coefs=coefs,
                control = co.x,
                constraint.tolerance = constraint.tolerance,
                print.level = print.level,
                base.weight = base.weight,
                full.t=full.t,
                full.c=full.c,
                shuffle.treat=shuffle.treat)

		if(!is.null(group.alpha)){
			min.c$solution[which(group.alpha==0)] <- 0
		}
		penalty <- rep(c(0, min.c$solution), times=grouping)
		if (!is.null(penalty.pos)){
			penalty[penalty.pos] <- penalty.val
		}

		z <- hb(
			tr_total=as.matrix(tr.total),
			co_x=co.x,
			coefs=as.matrix(coefs),
			base_weight=as.matrix(base.weight),
			alpha=penalty, 
			max_iterations=max.iterations,
			constraint_tolerance=constraint.tolerance,
			print_level=print.level
			)
	} #end of if, else


	weights<-z$Weights_ebal
	cc <- z$coefs

	if (ds){
		out.mat <- selected$resX
	}
	if (expand.degree>1 & !ds){
		out.mat <- expand$mat
	}
	if (expand.degree<=1 & !ds){
		out.mat <- X
	}

	grouping[1] <- grouping[1]-1
	group.assignment <- grouping

	out <- list(converged=z$converged,
				weights=weights, 
				coefs=cc, 
				Treatment=Treatment, 
				Y=unlist(data[,Y]), 
				mat=out.mat, 
				group.assignment=group.assignment,
				penalty=alpha,
				call=mcall)

	if (cv==TRUE){
		groups <- c("linear terms", "two-way interactions", "square terms", "three-way interactions", "linear-squre interactions", "cubic terms")
		out[["penalty"]] <- penalty[cumsum(out$group.assignment)+1]
		names(out[["penalty"]]) <- groups[1:length(out[["penalty"]])]
		out[["penalty"]][penalty.names] <- penalty.val
	}

	class(out) <- "hbal"

	return(out)
}
