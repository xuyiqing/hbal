#' @title Hierarchically Regularized Entropy Balancing
#' @aliases hbal
#' @description \code{hbal} performs hierarchically regularized entropy balancing 
#' such that the covariate distributions of the control group match those of the 
#' treatment group. \code{hbal} automatically expands the covariate space to include
#' higher order terms and uses cross-validation to select variable penalties for the 
#' balancing conditions.
#' @usage hbal(data, Treat, X, Y = NULL, w = NULL, 
#'      X.expand = NULL, X.keep = NULL, expand.degree = 1,
#'      coefs = NULL, max.iterations = 200, cv = NULL, folds = 4,
#'      ds = FALSE, group.exact = NULL, group.alpha = NULL,
#'      term.alpha = NULL, constraint.tolerance = 1e-3, print.level = 0,
#'      grouping = NULL, group.labs = NULL, linear.exact = TRUE, shuffle.treat = TRUE,
#'      exclude = NULL,force = FALSE, seed = 94035)
#' @param data                 a dataframe that contains the treatment, outcome, and covariates.   
#' @param Treat                a character string of the treatment variable.
#' @param X                    a character vector of covariate names to balance on.
#' @param Y                    a character string of the outcome variable.
#' @param w          		   a character string of the weighting variable for base weights
#' @param X.expand             a character vector of covariate names for serial expansion.
#' @param X.keep               a character vector of covariate names to keep regardless of whether they are selected in double selection.
#' @param expand.degree        degree of series expansion. 1 means no expansion. Default is 1.
#' @param coefs                initial coefficients for the reweighting algorithm (lambdas).
#' @param max.iterations       maximum number of iterations. Default is 200.
#' @param cv                   whether to use cross validation. Default is \code{TRUE}.
#' @param folds                number of folds for cross validation. Only used when cv is \code{TRUE}.
#' @param ds                   whether to perform double selection prior to balancing. Default is \code{FALSE}.
#' @param group.exact          binary indicator of whether each covariate group should be exact balanced.
#' @param group.alpha          penalty for each covariate group 
#' @param term.alpha           named vector of ridge penalties, only takes 0 or 1.
#' @param constraint.tolerance tolerance level for overall imbalance. Default is 1e-3.
#' @param print.level          details of printed output.
#' @param grouping             different groupings of the covariates. Must be specified if expand is \code{FALSE}.
#' @param group.labs           labels for user-supplied groups
#' @param linear.exact		   seek exact balance on the level terms
#' @param shuffle.treat        whether to use cross-validation on the treated units. Default is \code{TRUE}.
#' @param exclude              list of covariate name pairs or triplets to be excluded.
#' @param force                binary indicator of whether to expand covariates when there are too many
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
#' @importFrom stats setNames
#' @importFrom stats complete.cases
#' @importFrom stats weighted.mean
#' @importFrom stats density
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
	Y = NULL,
	w = NULL,
	X.expand = NULL, # variables to be expanded
	X.keep = NULL, # always keep despite double selection
	expand.degree = 1,
	coefs = NULL,
	max.iterations = 200,
	cv = NULL,
	folds = 4,
	ds = FALSE,
	group.exact = NULL,
	group.alpha = NULL,
	term.alpha = NULL,
	constraint.tolerance = 1e-3,
	print.level = 0,
	grouping = NULL,
	group.labs = NULL,
	linear.exact = TRUE,
	shuffle.treat = TRUE,
	exclude = NULL,
	force = FALSE,
	seed = 94035
	){

	# ntreated: number of treated units
	# ncontrols: number of control units # nolint # nolint
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

	
	#renames the covariates
	new_names <- paste0("X", seq_along(X))
	colnames(data)[colnames(data) %in% X] <- new_names
	if (is.null(X.expand) == FALSE)
	{
	  mapping <- setNames(new_names, X)
	  mapped_list <- mapping[match(X.expand, names(mapping))]
	  X.expand <- unname(mapped_list)
	}
	if (is.null(X.keep) == FALSE)
	{
	  mapping <- setNames(new_names, X)
	  mapped_list <- mapping[match(X.keep, names(mapping))]
	  X.keep <- unname(mapped_list)
	}
	old_names = X
	X = new_names
	
	# covariates for expansion
	if (expand.degree > 1 && is.null(X.expand) == TRUE) {
		if (print.level >= 1) message("All variables will be serially expanded")
		if ((expand.degree == 2 & length(X)>30) | (expand.degree == 3 & length(X)>10)) {
			if (force == TRUE) {
				message("Too many variables to expand. This can take up significant memory.\nPlease consider including fewer covariates or using the \`X.expand\` option\n")
			} else {
				stop("Too many variables to expand. This can take up significant memory.\nPlease consider including fewer covariates or using the \`X.expand\` option.\nAdd \'force = TRUE\' if you want to keep going")
			}
		}
		X.expand <- X
	}
	if (expand.degree == 1 && is.null(X.expand) == FALSE) {
		if (print.level >= 1) {
			message("\"X.expand\" is ignored because expand.degree = 1; no serial expansion will be performed")}
		X.expand <- NULL
	}

	# combine X, X.keep, X.levelonly
	X.all <- unique(c(X, X.keep, X.expand))
	X.levelonly <- setdiff(X.all, X.expand)
	
	# check if variables are in the dataset
	var.names <- colnames(data)
	if(!all(c(Treat, X.all, Y, w) %in% var.names)) stop("Some variable(s) specified are not in the data")
	
	valid_cols <- !sapply(data[, c(Treat, X.all, Y, w)], class) == 'character' # need all numeric columns
	if (sum(valid_cols) != length(valid_cols)) {
		stop('Some columns in the data are character columns. Consider converting them to numeric')
	}

	# listwise deletion
	valid_rows <- complete.cases(data[, c(Treat, X.all, Y, w)])
	if (sum(valid_rows) < nrow(data)) {
		if (print.level >= 1) {message("Some rows are dropped because they contain missing/NA/infinite values\n")}
		data <- data[valid_rows,]		
	}

	# Treatment indicator
	Treatment <- unname(unlist(data[,Treat]))
	names(Treatment) <- Treat
	if (sd(Treatment) == 0) stop("No variation in the treatment variable.")
	if (!is.numeric(Treatment)) stop("Treatment variable needs to be numeric")
	ntreated  <- sum(Treatment==1)
	ncontrols <- sum(Treatment==0)  

	# base weights
	if (is.null(w)==FALSE) {
		base.weights <- data[, w]
	} else {
		base.weights <- rep(1, nrow(data))
	}
	base.weights.tr <- base.weights[Treatment==1]
	base.weights.co <- base.weights[Treatment==0]
		
	###########################
	# Covariates X
	###########################

	## logic ##
	# 1. scale
	# 2. check collinearity, drop
	# 3. serial expansion on original data, then scale
	# 4. check collinearity again, drop
	# 5. double selection

	# Data (Matrices): X, X.sav
	# Variable names: X.all, X.keep, X.levelonly
	# Variable positions: X.levelonly.pos, X.expand.pos, X.keep.pos

	# check X variation
	X.novar <- X.all[which(apply(data[, X.all, drop = FALSE], 2, sd) == 0)] # controls of no variations
	if (length(X.novar) > 0 && print.level >= 0) {
		message(paste("The following variable(s) have no variations and are automatically dropped:", paste0(X.novar, collapse = ", "),"\n"))
		X.all <- setdiff(X.all, X.novar)
	}

	# save a copy of data (to be saved in the output)
	X.sav <- as.matrix(data[, X.all, drop = FALSE])

	# scaling
	X <- scale(X.sav)

	# check collinearity
	if (qr(X)$rank != ncol(X)) {
		message("Collinearity in covariate matrix for controls; ")
		X.good.pos <- 1
		for (i in 2:ncol(X)) {
			X.test.pos <- c(X.good.pos, i)
			if (qr(X[,X.test.pos])$rank == length(X.test.pos)) {
				X.good.pos <- c(X.good.pos, i)
			}
		}
		X.bad.pos <- setdiff(1:ncol(X), X.good.pos)
		X.bad <- X.all[X.bad.pos]
		if (print.level >= 0) message("the following variable(s) are removed:",paste(X.bad, collapse = ", "),"\n")
		# update X & saved X data (not scaled)
		X <- X[, X.good.pos, drop = FALSE]
		X.sav <- X.sav[, X.good.pos, drop = FALSE]
		# update variable names
		X.all <- X.all[X.good.pos]
		X.keep <- setdiff(X.keep, X.bad)		
		X.levelonly <- setdiff(X.levelonly, X.bad)		
	}
	
	# get column numbers
	X.levelonly.pos <- which(X.all %in% X.levelonly)
	X.expand.pos <- which(X.all %in% X.expand)
	X.keep.pos <- which(X.all %in% X.keep)

	# grouping and covariate expansion	
	if(is.null(grouping)==FALSE) {
		if (expand.degree >1) {stop("User-supplied groups using \"grouping\" are incompatible with serial expansion")}
		if (length(grouping)!=length(group.labs)) {stop(paste("Please supply",length(grouping),"group labels"))}
		names(grouping) <- group.labs
	} else {
		grouping <- ncol(X)
		names(grouping) <- "linear"
	}
	if(length(expand.degree)!=1) {stop("\"expand.degree\" needs to be of length 1")}
	if (!expand.degree %in% c(1, 2, 3)) {stop("\"expand.degree\" needs to be one of c(1, 2, 3)")}

	# series expansion of the covariates
	if (expand.degree > 1) {
		expand <- covarExpand(X.sav[, X.expand.pos, drop = FALSE], 
			exp.degree = expand.degree, treatment = Treatment, exclude = exclude) 
		X.tmp <- expand$mat # expanded covariates in matrix form
		grouping <- expand$grouping # grouping of the expanded covariates
		if (expand.degree == 2) {group.labs <- c("linear", "squared", "two-way")} # three of them
		if (expand.degree == 3) {group.labs <- c("linear", "two-way", "squared", "three-way", "squared*linear", "cubic")} # six of them
		names(grouping) <- group.labs[1:length(grouping)]
		## add level-only covariates
		if (length(X.levelonly.pos)>0) { 
			X.sav <- cbind(X.sav, X.tmp[, -c(1:length(X.expand.pos))])
			grouping[1] <- grouping[1] + length(X.levelonly.pos)
		} else {
			X.sav <- X.tmp
		}
		# remove empty
		if (0 %in% grouping){
			grouping <- grouping[-which(grouping==0)]
		}

		# X.sav is unscaled; X is scaled
		X <- scale(X.sav)

		# check collinearity again after serial expansion
		grouping.expand <- c()
		for (i in 1:length(grouping)) {
			grouping.expand <- c(grouping.expand, rep(names(grouping)[i], grouping[i]))
		}
		if (qr(X)$rank != ncol(X)) {
			X.good.pos <- 1
			for (i in 2:ncol(X)) {
				X.test.pos <- c(X.good.pos, i)
				if (qr(X[,X.test.pos])$rank == length(X.test.pos)) {
					X.good.pos <- c(X.good.pos, i)
				}
			}
			X.bad.pos <- setdiff(1:ncol(X), X.good.pos) # these are column numbers
			if (print.level >= 0) {
				message("After serial expansion, the following variable(s) are removed due to collinearity:",
					paste(colnames(X)[X.bad.pos], collapse = ", "),"\n")}
			X <- X[, X.good.pos, drop = FALSE]
			grouping.expand <- grouping.expand[-X.bad.pos]
			grouping <- as.numeric(table(grouping.expand)) # with grouping names
			names(grouping) <- names(table(grouping.expand))
			# update X.sav
			X.sav <- X.sav[, X.good.pos, drop = FALSE]
			# update X.keep.pos
			X.keep.intersect <- intersect(X.keep.pos, X.good.pos)
			X.keep.pos <- which(X.good.pos %in% X.keep.intersect)
		}
	} 

	full <- X # need to keep a copy of the standardized X for cross-validation. Normalization so that means and variances of covariates on the same scale 
	full.t <- full[Treatment==1,]
	full.c <- full[Treatment==0,]

	# double selection
    if (ds == TRUE){
		selected <- doubleSelection(X=X, W=Treatment, Y=unlist(data[,Y]), 
			grouping=grouping)
		X.pos <- X.select.pos <- selected$covar.keep
		grouping <- selected$penalty.list	
		# add those must keep	
		if (length(X.keep.pos)>0) {
			X.pos <- union(X.keep.pos, X.select.pos)
			grouping[1] <- grouping[1] + (length(X.pos) - length(X.select.pos))			
		}
		# update X and X.sav
		X <- X[, X.pos, drop = FALSE]		
		X.sav <- X.sav[, X.pos, drop = FALSE] 		
	}

	full.c <- cbind(rep(1,ncontrols), full.c)
	co.x <- X[Treatment==0,,drop = FALSE] # control group
	co.x <- cbind(rep(1,ncontrols),co.x)
	
	# treated group mean
	tr.total <- c(1, apply(X[Treatment==1,,drop=FALSE], 2, weighted.mean, w = base.weights.tr))	
	grouping[1] <- grouping[1] + 1 # need one more for normalizing constraint


	###########################
	# User-Supplied Penalties
	###########################

	# penalty on specific terms
	if (!is.null(term.alpha)){
	  if (expand.degree > 1) {
	    warning("\"term.alpha\ does not work with feature expansion; ignored.")
	    term.alpha <- NULL
	  }
	}
	if (!is.null(term.alpha)){
	 	penalty.names <- names(term.alpha)
		if (length(match(penalty.names, colnames(X)))==0) {
			if (print.level >= 1) message("Invalid variable name(s); \"term.alpha\" is ignored\n")
			penalty.names <- penalty.val <- penalty.pos <- NULL
		} else {
			penalty.val <- term.alpha
			penalty.pos <- match(penalty.names, old_names)+1
		}		
	}else{
		penalty.names <- penalty.val <- penalty.pos <- NULL
	}

	##############
	# Checks
	##############
	
	if (sum(grouping) != ncol(co.x)) stop("Sum of grouping must be equal to ncol(X)")
	if (sum(Treatment != 1 & Treatment != 0) > 0) stop("Treatment indicator ('Treatment') must be a logical variable, TRUE (1) or FALSE (0)")
	if (var(Treatment) == 0) stop("Treatment indicator ('Treatment') must contain both treatment and control observations")
	if (sum(is.na(Treatment)) > 0) stop("Treatment contains missing data")
	if (length(Treatment) != nrow(X)) stop("length(Treatment) != nrow(X)")
	if (sum(is.na(X)) > 0) stop("X contains missing data")
	if (length(max.iterations) != 1 ) stop("length(max.iterations) != 1")
	if (length(constraint.tolerance) != 1 ) stop("length(constraint.tolerance) != 1")
	if (qr(co.x)$rank != ncol(co.x)) stop("collinearity in covariate matrix for controls (remove collinear covariates)")
	
	# initialization for lambdas
	if (is.null(coefs)) {coefs = c(log(tr.total[1]/sum(base.weights.co)),rep(0,(ncol(co.x)-1)))}
	if (length(coefs) != ncol(co.x)) {stop("coefs needs to have same length as number of covariates plus one")}

	# Prepare for cross-validation
	if (is.null(cv) == FALSE) {
		if (cv %in% c(TRUE, FALSE) == FALSE) {stop("cv needs to be of type logical")}
		# if (cv==TRUE && length(grouping)==1){
		# 	cv <- FALSE
		# 	if (print.level >= 1) {
		# 		message("length(grouping)==1, no need to cross-validate tuning parameters. \n
		# 			Either double selection selected 0 higher order term or the supplied grouping has length 1")}
		# }
	} else {
		cv <- FALSE
	}
	
	
	# check group.alpha; if correct, no need to cv
	if (is.null(group.alpha)==FALSE) {
		if (length(group.alpha) != length(grouping)) {
			stop(paste("\"group.alpha\" should be of the same length as the number of groups, which is",length(grouping)))
		} 
		if (sum(group.alpha<0)>0) {
			stop("Elements in \"group.alpha\" should be non-negative")
		}
		if (cv == TRUE) {
			if (print.level >= 0) message("Cross-validation is skipped when \"group.alpha\" is supplied\n")
			cv <- FALSE
		} 
	} else {
		if (cv == FALSE) {
			group.alpha <- rep(0, length(grouping)) # no group.alpha, no cv
		}
	}


	# check group.exact
	if (is.null(group.exact)==FALSE) {
		if (cv == FALSE) {
			group.exact <- NULL
			if (print.level >= 0) {
				message("\"group.exact\" is ignored when \"cv = FALSE\" because either all covariates 
					will be exactly balanced on or \"group.alpha\" is supplied\n")}
		} else {
			if (sum(group.exact %in% c(0,1))!= length(group.exact)) {
				stop("Elements in \"group.exact\" should be TRUE or FALSE")
			}
		}
	}	

	if (linear.exact == TRUE) {
		if (is.null(group.exact) == TRUE) {
			group.exact <- c(1, rep(0, length(grouping)-1))
		} else {
			group.exact[1] <- 1
		}
	}


	# Cross-validation
	if (cv==TRUE){
		message("Crossvalidation...\n")
		fold.num.co <- rep(1:folds, ceiling(ncontrols/folds))
		fold.co <- sample(fold.num.co, ncontrols, replace=FALSE) 
		fold.num.tr <- rep(1:folds, ceiling(ntreated/folds))
		fold.tr <- sample(fold.num.tr, ntreated, replace=FALSE)

		min.c <- nloptr(x0 = rep(0, length(grouping)),
			eval_f = crossValidate,
			lb = rep(0, length(grouping)),
			ub = rep(100, length(grouping)),
			opts = list('algorithm'='NLOPT_LN_COBYLA',
						'maxeval' =200,
			'ftol_rel'=1e-3,
			'ftol_abs'=1e-5,
			'xtol_abs'=1e-3,
			'print_level'=print.level),
			penalty.pos=penalty.pos,
			penalty.val=penalty.val,
			group.exact=group.exact,
			grouping=grouping,
			folds=folds,
			treatment = X[Treatment==1,],
			fold.co = fold.co,
			fold.tr=fold.tr,
			coefs=coefs,
			control = co.x,
			constraint.tolerance = constraint.tolerance,
			  print.level = print.level,
			base.weight = base.weights.co,
			full.t=full.t,
			full.c=full.c,
			shuffle.treat=shuffle.treat)

		names(min.c)		
		min.c$lower_bounds
		min.c$iterations
		min.c$message
		min.c$solution		
		min.c$iterations
		min.c$objective

		if(!is.null(group.exact)){
			min.c$solution[which(group.exact==1)] <- 0
		}
		# save CV results: penalty (same length of "grouping")
		group.alpha <- min.c$solution		
		
	} #end of tuning

	# expand to penalty for individual terms
	names(group.alpha) <- names(grouping)
	alpha <- rep(group.alpha, times=grouping) 
	# apply individual penalty values for specific terms
	if (!is.null(penalty.pos)) { # fill in specific alpha for individual terms
		alpha[penalty.pos] <- penalty.val
	}
	alpha[1]  <- 0 # normalizing term should not have penalty		

	##################
	# Main Algorithm
	##################

	z <- hb(
			tr_total=as.matrix(tr.total),
			co_x=co.x,
			coefs=as.matrix(coefs),
			base_weight=as.matrix(base.weights.co),
			alpha=alpha, 
			max_iterations=max.iterations,
			constraint_tolerance=constraint.tolerance,
			print_level=print.level
			)
  
	# non converge warning
	if(!z$converged)
	{
	  message("Not converged. Try setting ds = TRUE or using X.keep and exclude option.\n")
	}
	
	# save weights and coefficients
	

	weights <- rep(NA, length(Treatment))
	weights.co <- z$Weights_ebal
	weights.co <- weights.co/sum(weights.co)*sum(base.weights.tr) # normalize to the total number of treated
	weights[Treatment == 0] <- weights.co
	weights[Treatment == 1] <- base.weights.tr
	cc <- z$coefs

	# take out the normalizing constant
	grouping[1] <- grouping[1]-1 
	alpha <- alpha[-1]
	Covar <- colnames(X.sav)
	#rename the covar
	for(i in 1:length(new_names)) Covar <- stringr::str_replace_all(Covar, paste0("\\b",new_names[i], "\\b"), str_trunc(old_names[i], 5, side="right", ellipsis=""))
	Covar_nonum <- Covar
	Covar <- paste(Covar, seq_along(Covar), sep = ".")
	colnames(X.sav) <- Covar
	names(alpha) <- Covar
	
	## balance table
	bal.tab <- matrix(NA, length(Covar), 5) # tr, co, co.w, diff, diff.w
	rownames(bal.tab) <- Covar_nonum
	treat <- X.sav[Treatment==1, , drop = FALSE] # treated
	control <- X.sav[Treatment==0, , drop = FALSE] # control
	bal.tab[,1] <- apply(treat, 2, weighted.mean, w = base.weights.tr)
	bal.tab[,2] <- apply(control, 2, weighted.mean, w = base.weights.co)	 
	bal.tab[,3] <- apply(control, 2, weighted.mean, w = weights.co)	 
	denom <- apply(X.sav, 2, sd)
	bal.tab[,4] <- (bal.tab[,1] - bal.tab[,2])/denom
	bal.tab[,5] <- (bal.tab[,1] - bal.tab[,3])/denom
	colnames(bal.tab) <- c("Tr.Mean", "Co.Mean", "W.Co.Mean", 
	"Std.Diff.(O)", "Std.Diff.(W)")
	bal.tab <- round(bal.tab, 2)
	
	out <- list(converged=z$converged,
				weights=weights, 
				weights.co=weights.co,
				coefs=cc, 
				Treatment = Treatment,				
				mat=X.sav, # expanded covariates
				grouping=grouping,
				group.penalty=group.alpha,
				term.penalty=alpha,
				bal.tab = bal.tab,
				base.weights = base.weights,
				Treat = Treat)	

	if (is.null(Y)==FALSE) {
		out <- c(out, list(Outcome = data[,Y], Y = Y))
	}
	out <- c(out, call = mcall)

	class(out) <- "hbal"

	return(out)
}
