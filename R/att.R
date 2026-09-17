#' @title Estimating the ATT from an hbal object
#' @aliases att
#' @description \code{att} estimates the average treatment effect on the treated (ATT) from an 
#' hbal object returned by \code{hbal}. 
#' @usage att(hbalobject, method="abw", dr=TRUE, displayAll=FALSE, alpha=0.9,
#'     seed=NULL, nfolds=5, ...)
#' @param hbalobject  an object of class \code{hbal} as returned by \code{hbal}.
#' @param method      estimation method for the ATT. The default \code{"abw"}
#'   (Augmented Balancing Weights; new in version 1.3.0) is a cross-fitted,
#'   Neyman-orthogonal estimator that combines the hbal weights with a
#'   ridge-regularized control-group outcome model (see Details).
#'   \code{"lm_robust"} reproduces the default of versions before 1.3.0 (weighted
#'   regression via \code{lm_robust}); \code{"lm_lin"} (Lin 2013) and
#'   \code{"elnet"} (Athey, Imbens and Wager 2018) are also available.
#' @param dr          doubly robust, whether an outcome model is included in estimating
#'   the ATT. With \code{dr = FALSE} every method returns the weighted difference in
#'   means. \code{"abw"} and \code{"lm_robust"} compute it directly; \code{"lm_lin"} and
#'   \code{"elnet"} are defined only with an outcome model, so \code{att} estimates them
#'   as \code{"lm_robust"}, giving the same estimate, standard error and degrees of
#'   freedom as \code{att(x, method = "lm_robust", dr = FALSE)}.
#' @param displayAll  only displays treatment effect by default. If \code{TRUE},
#'   returns the fitted \code{lm_robust} or \code{lm_lin} object, or, for
#'   \code{method = "abw"}, a list with the estimate and its components (see Value).
#' @param alpha       tuning parameter for glmnet (\code{method = "elnet"} only).
#' @param seed        a single number that fixes the cross-fitting fold assignment of
#'   \code{method = "abw"}; \code{NULL} (the default) assigns folds in the row order
#'   of the controls. The assignment is a deterministic function of \code{seed};
#'   \code{set.seed()} is never called. Ignored by the other methods.
#' @param nfolds      number of cross-fitting folds for \code{method = "abw"}. Default
#'   is 5. Reduced automatically, with a message, when the control group is small
#'   relative to the number of covariates (see Details). Ignored by the other methods.
#' @param ...         arguments passed to lm_lin or lm_robust (e.g. \code{se_type},
#'   \code{clusters}). Not accepted by \code{method = "abw"}.
#' @details The default \code{method = "abw"} (Augmented Balancing Weights;
#'   Ben-Michael, Feller, Hirshberg and Zubizarreta 2021; Bruns-Smith, Dukes, Feller
#'   and Ogburn 2023) combines the hbal weights with a control-group outcome model in
#'   an augmented moment condition that is Neyman-orthogonal: its first-order
#'   sensitivity to estimation error in either nuisance component (the outcome model,
#'   given hbal's balance on \code{mat}; the weights, given a correct outcome model)
#'   is zero. Let \eqn{w_i} be the base weights of the treated units, \eqn{W_1}{W1}
#'   their sum, \eqn{\gamma_i}{gamma_i} the hbal weights of the controls
#'   (\code{weights.co}, which also sum to \eqn{W_1}{W1}), and \eqn{\hat{\mu}_0}{mu0}
#'   a ridge regression of the outcome on the columns of \code{hbalobject$mat},
#'   fitted on the controls only and weighted by their \eqn{\gamma_i}{gamma_i}:
#'   every column is standardized with the mean and
#'   standard deviation of the controls used in the fit, the intercept is not
#'   penalized, and the penalty \eqn{\lambda}{lambda} is chosen once per call from
#'   the full control sample by generalized cross-validation (GCV; Golub, Heath and
#'   Wahba 1979) over a fixed grid of candidate values, restricted to candidates
#'   whose effective degrees of freedom do not exceed
#'   \eqn{\min(p, n_0 / 2) + 1}{min(p, n0 / 2) + 1}, with \eqn{n_0}{n0} controls
#'   and \eqn{p} columns in \code{mat}, so that the fit can never approach an
#'   interpolation of the controls. The selection is deterministic: no data
#'   splitting and no random numbers are involved, and \eqn{\lambda = 0}{lambda = 0}
#'   (weighted least squares) is chosen whenever it is admissible and minimizes the
#'   GCV score. The estimate is
#'   \deqn{\hat{\tau} = \frac{1}{W_1}\Big[\sum_{T_i = 1} w_i (Y_i - \hat{m}_i) -
#'   \sum_{T_i = 0} \gamma_i \hat{e}_i\Big],}{tau = (1 / W1) [ sum_{T=1} w_i (Y_i - m_i) -
#'   sum_{T=0} gamma_i e_i ],}
#'   where \eqn{\hat{m}_i}{m_i} is the outcome model fitted on all controls and evaluated
#'   at treated unit \eqn{i}, and \eqn{\hat{e}_i}{e_i} is the cross-fitted residual of
#'   control \eqn{i}: the controls are split into \code{nfolds} folds and each control's
#'   residual uses the fit obtained without its own fold. The number of folds actually
#'   used is \eqn{K = \max(1, \min(\mathrm{nfolds}, \lfloor n_0 / (2p) \rfloor))}{K = max(1, min(nfolds, floor(n0 / (2p))))},
#'   (\eqn{n_0}{n0} and \eqn{p} as above); \eqn{K = 1} means no
#'   cross-fitting. Setting \code{nfolds = 1} disables cross-fitting (the
#'   correction term above becomes identically zero) and lowers root-mean-squared
#'   error at small sample sizes in simulations accompanying this release, at the
#'   cost of an understated standard error and, under outcome-model
#'   misspecification, added bias; \code{nfolds = 5} (the default) is recommended
#'   unless a smaller root-mean-squared error is wanted more than a reliable
#'   standard error. The standard error is influence-function based: with
#'   \eqn{\psi_i = w_i T_i (Y_i - \hat{m}_i - \hat{\tau}) - \gamma_i (1 - T_i) \hat{e}_i}{psi_i = w_i T_i (Y_i - m_i - tau) - gamma_i (1 - T_i) e_i},
#'   the variance estimate is \eqn{\sum_i \psi_i^2 / W_1^2}{sum_i psi_i^2 / W1^2}, treating
#'   the weights and the outcome fit as fixed; p-values and the 95 percent confidence
#'   interval use a t distribution with \eqn{n - 1} degrees of freedom. With
#'   \code{dr = FALSE} the outcome model is dropped (\eqn{\hat{m}_i = 0}{m_i = 0},
#'   \eqn{\hat{e}_i = Y_i}{e_i = Y_i}) and the same formulas give the weighted difference
#'   in means. Fold assignment is deterministic: repeated calls with the same
#'   \code{seed} (including \code{NULL}) return identical results. If the
#'   entropy-balancing weights of \code{hbalobject} did not converge
#'   (\code{hbalobject$converged} is \code{FALSE}), \code{att} issues a warning for
#'   either value of \code{dr}, because the estimate and its standard error may then
#'   be unreliable.
#'
#'   \code{method = "lm_robust"} (the default before version 1.3.0) and
#'   \code{method = "lm_lin"} are wrappers for \code{lm_robust} and \code{lm_lin} from
#'   the \pkg{estimatr} package, fitted with the hbal weights; \code{method = "elnet"}
#'   implements the approximate residual balancing estimator of Athey, Imbens and Wager
#'   (2018) via \pkg{glmnet}. For \code{method = "lm_robust"}, the default
#'   \code{se_type} understates the standard error by 5 to 7 percent in the same
#'   simulations; \code{se_type = "HC3"} corrects this for samples of 500 or fewer
#'   but still runs about 5 percent short at 1,000.
#'   Before version 1.3.0, \code{dr = FALSE} was accepted and silently ignored by
#'   \code{method = "lm_lin"} and \code{method = "elnet"}; both now fall back to
#'   \code{"lm_robust"}.
#' @references Ben-Michael, E., Feller, A., Hirshberg, D. A., and Zubizarreta, J. R.
#'   (2021). The balancing act in causal inference. arXiv:2110.14831.
#'
#'   Bruns-Smith, D., Dukes, O., Feller, A., and Ogburn, E. L. (2023). Augmented
#'   balancing weights as linear regression. arXiv:2304.14545.
#'
#'   Golub, G. H., Heath, M., and Wahba, G. (1979). Generalized cross-validation
#'   as a method for choosing a good ridge parameter. Technometrics, 21(2),
#'   215-223.
#'
#'   Athey, S., Imbens, G. W., and Wager, S. (2018). Approximate residual
#'   balancing: debiased inference of average treatment effects in high
#'   dimensions. Journal of the Royal Statistical Society: Series B, 80(4),
#'   597-623.
#'
#'   Lin, W. (2013). Agnostic notes on regression adjustments to experimental
#'   data: Reexamining Freedman's critique. The Annals of Applied Statistics,
#'   7(1), 295-318.
#' @return A data frame with one row and seven columns (Estimate, Std. Error,
#'   t value, Pr(>|t|), CI Lower, CI Upper, DF) when \code{displayAll = FALSE}
#'   (the default) or \code{method = "elnet"}. When \code{displayAll = TRUE}: for
#'   \code{method = "lm_robust"} or \code{"lm_lin"}, the fitted model object returned
#'   by \code{lm_robust} or \code{lm_lin}; for \code{method = "abw"}, a plain list
#'   with elements
#'   \describe{
#'     \item{estimate}{the ATT point estimate.}
#'     \item{se}{its influence-function standard error.}
#'     \item{df}{degrees of freedom of the t distribution used for the p-value and
#'     the confidence interval, \eqn{n - 1}.}
#'     \item{method}{\code{"abw"}.}
#'     \item{dr}{the \code{dr} value used.}
#'     \item{influence}{numeric vector of length \eqn{n}: the per-unit moment
#'     contributions \eqn{\psi_i}{psi_i} (see Details), in the row order of
#'     \code{hbalobject}; they sum to zero.}
#'     \item{fold}{integer vector of length \eqn{n}: the cross-fitting fold of each
#'     control (\code{NA} for treated units, and for every unit when \code{dr = FALSE}).}
#'     \item{nfolds}{the number of folds actually used (\code{NA} when \code{dr = FALSE}).}
#'     \item{nuisance}{a list with \code{coef_full} (coefficients of the ridge outcome
#'     model fitted on all controls, on the original scale of \code{mat}, named
#'     \code{"(Intercept)"} followed by the columns of \code{mat}), \code{coef_folds}
#'     (a list of the \code{nfolds} fold-specific coefficient vectors), and
#'     \code{rank_full} (the rank of the full-control design including the
#'     intercept); \code{NULL}, \code{list()} and \code{NA} respectively when
#'     \code{dr = FALSE}.}
#'     \item{weights}{a list with \code{treated} (the base weights of the treated
#'     units) and \code{control} (the hbal weights of the controls), each in the row
#'     order of \code{hbalobject$mat} within its group.}
#'   }
#' @importFrom estimatr lm_lin lm_robust
#' @importFrom stats as.formula
#' @importFrom stats predict
#' @importFrom stats pt
#' @importFrom stats qt
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
#' sout <- summary(att(out))       # default: augmented balancing weights (abw)
#' att(out, method = "lm_robust")  # the default before version 1.3.0
#' @export

# Internal helper for att(): given whatever generics::tidy() returned on an
# estimatr fit (a plain data.frame under estimatr < 2.0.0, a tibble under
# estimatr >= 2.0.0), select the treatment row and the seven display columns
# by name, and return a plain data.frame with the package's display column
# names and the treatment name as its row name. Not exported.
# @noRd
.att_tidy_select <- function(tidy_out, Tr) {
	tidy_out <- as.data.frame(tidy_out, stringsAsFactors = FALSE)

	row_idx <- which(tidy_out[["term"]] == Tr)
	if (length(row_idx) != 1) {
		stop(
			"att(): expected exactly one row with term == \"", Tr, "\" in the ",
			"tidy() output from estimatr, found ", length(row_idx), ". ",
			"This usually means the treatment term name changed or estimatr's ",
			"tidy() output is not in the expected shape."
		)
	}

	required_cols <- c("estimate", "std.error", "statistic", "p.value",
						"conf.low", "conf.high", "df")
	missing_cols <- setdiff(required_cols, colnames(tidy_out))
	if (length(missing_cols) > 0) {
		stop(
			"att(): the tidy() output from estimatr is missing expected ",
			"column(s): ", paste(missing_cols, collapse = ", "), ". ",
			"This usually means an incompatible estimatr version is installed."
		)
	}

	out <- tidy_out[row_idx, required_cols, drop = FALSE]
	colnames(out) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)",
						"CI Lower", "CI Upper", "DF")
	rownames(out) <- Tr
	out
}

att <- function(
	hbalobject,
	method="abw",
	dr=TRUE,
	displayAll=FALSE,
	alpha = 0.9,
	seed = NULL,
	nfolds = 5,
	...
	){
	if(!inherits(hbalobject, "hbal")){
		stop("hbalobject must be an hbal object from a call to hbal()")
    }
    elpss <- list(...)
	# "lm_lin" and "elnet" are defined only with an outcome model, so under
	# dr = FALSE they fall back to "lm_robust", whose dr = FALSE branch is the
	# weighted difference in means. "abw" and "lm_robust" consult dr themselves
	# and must not be rerouted. The old test compared method against "lin", a
	# string no branch below uses, so the fallback never fired (GitHub PR #8).
	if (dr == FALSE && method %in% c("lm_lin", "elnet")) {
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
	# Cross-fitted, Neyman-orthogonal augmented estimator (default since 1.3.0).
	# Self-contained: builds its own output in R/att_abw.R and never reaches the
	# estimatr / tidy() tail below. seed and nfolds are validated here only; the
	# legacy methods ignore them.
	if (method == "abw") {
		if (length(elpss) > 0) {
			dots_names <- names(elpss)
			if (is.null(dots_names)) dots_names <- rep("", length(elpss))
			dots_names[dots_names == ""] <- "<unnamed>"
			stop(sprintf("att(): method = \"abw\" does not accept: %s. These arguments (e.g. se_type, clusters) apply only to method = \"lm_robust\" or \"lm_lin\".",
				paste(dots_names, collapse = ", ")))
		}
		if (!(length(nfolds) == 1 && is.numeric(nfolds) && is.finite(nfolds) && nfolds == round(nfolds) && nfolds >= 1)) {
			stop("att(): \"nfolds\" must be a single positive whole number")
		}
		if (!is.null(seed) && !(length(seed) == 1 && is.numeric(seed) && is.finite(seed))) {
			stop("att(): \"seed\" must be NULL or a single finite number")
		}
		res <- .abw_att(hbalobject, dr = dr, seed = seed, nfolds = nfolds)
		if (displayAll == FALSE) {
			return(.abw_table(res, Tr))
		} else {
			return(res)
		}
	}
	# linear regression
	else if (method=="lm_robust"){
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
	else if (method=="lm_lin"){
		covariates <- as.formula(paste0("~ ", paste0(Covar, collapse=" + ")))
		if (is.null(elpss[['se_type']])){
			out <- lm_lin(formula = as.formula(paste0(Y, ' ~ ', Tr)), covariates=covariates, weights=w, data=dat, se_type='stata', ...)
		}else{
			out <- lm_lin(formula = as.formula(paste0(Y, ' ~ ', Tr)), covariates=covariates, weights=w, data=dat, ...)
		}
	}
	# Athey (2018): Approximate residual balancing
	else if (method == "elnet"){
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
	else {
		stop(sprintf("att(): unrecognized method \"%s\". Valid values are \"abw\", \"lm_robust\", \"lm_lin\", \"elnet\".", method))
	}
	# DisplayAll Option of ATT
	if (displayAll == FALSE){
		out <- .att_tidy_select(tidy(out), Tr)
		return(out)
	}
	else
	{
	  return(out)
	}
	
}