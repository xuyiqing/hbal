# Internal helpers for att(method = "abw"): the cross-fitted, Neyman-orthogonal
# augmented balancing weights (ABW) ATT estimator that is the default of att()
# since hbal 1.3.0. Nothing in this file is exported, and no function here
# reads or writes .Random.seed (no set.seed(), sample(), runif(), ...): the fold
# assignment is a deterministic hash of the `seed` argument, and the ridge
# penalty of the outcome model is a deterministic plug-in (no cross-validation).
#
# Notation (i indexes the rows of the hbal object in their original order):
#   T_i      treatment indicator                        hbalobject$Treatment
#   X_i      covariate row, unscaled, already expanded   hbalobject$mat
#   Y_i      outcome                                     hbalobject$Outcome
#   w_i      weight of every unit                        hbalobject$weights
#            (user-supplied w or 1 for the treated; the balancing weight for
#            the controls)
#   gamma_i  balancing weight of control i               hbalobject$weights.co
#            (sums to W1 = sum of the treated base weights)
#   m_i      = mu0_full(X_i): the control-only outcome model, fitted on ALL
#            controls, evaluated at treated unit i
#   e_i      = Y_i - mu0_(-k(i))(X_i): the cross-fitted residual of control i,
#            mu0_(-k) being fitted on the controls outside i's fold k(i)
#
# Point estimate
#   tau = (1 / W1) * [ sum_{T_i=1} w_i (Y_i - m_i) - sum_{T_i=0} gamma_i e_i ]
# Per-unit moment contribution (sums to zero at tau by construction)
#   psi_i = w_i T_i (Y_i - m_i - tau) - gamma_i (1 - T_i) e_i
# Variance (the bread d/dtau sum_i psi_i = -W1 is a known constant)
#   V = sum_i psi_i^2 / W1^2
# dr = FALSE sets m_i = 0 and e_i = Y_i in the same formulas, which reduces
# tau to the weighted difference in means with the matching sandwich variance.
#
# Outcome model (every fit, full sample and each fold): weighted ridge
# regression of Y on the columns of X, each column standardized with the
# training rows' own mean and standard deviation, intercept unpenalized:
#   beta_s = argmin sum_i w_i (Y_i - [1, x_s_i] beta)^2 + lambda * ||beta[-1]||^2
# lambda >= 0 is computed ONCE per att() call from the full control sample by
# the Hoerl-Kennard (1970) plug-in lambda = p * sigma2 / ||beta_ols[-1]||^2
# (sigma2 = weighted residual variance of the lambda = 0 fit) and reused for
# every fit. lambda = 0 reproduces weighted least squares exactly.

# Multiplicative-hash constants (Knuth / xxHash primes) used to turn `seed`
# into a deterministic permutation of the controls without touching R's RNG.
.ABW_HASH_A <- 2654435761
.ABW_HASH_B <- 2246822519
.ABW_HASH_MOD <- 4294967296  # 2^32
# Upper quantile of the t distribution for the two-sided 95 percent interval.
.ABW_CI_PROB <- 0.975
# Cross-fitting keeps at least this many controls per fold per covariate
# column, so that fold-level fits stay reasonably determined:
# K <= floor(n0 / (.ABW_CONTROLS_PER_FOLD_PER_COVARIATE * p)).
.ABW_CONTROLS_PER_FOLD_PER_COVARIATE <- 2
# Defensive cap on the Hoerl-Kennard ridge penalty (standardized design), and
# the threshold below which the OLS slope norm counts as degenerate (lambda
# then falls back to 0, i.e. plain weighted least squares).
.ABW_LAMBDA_MAX <- 50
.ABW_BETA_NORM_TOL <- 1e-8

#' Deterministic permutation of the controls used for fold assignment
#'
#' @param n0 number of controls.
#' @param seed \code{NULL} (identity order) or a single finite number.
#' @return an integer permutation of \code{seq_len(n0)}.
#' @keywords internal
#' @noRd
.abw_hash_order <- function(n0, seed) {
	j <- as.numeric(seq_len(n0))
	if (is.null(seed)) {
		h <- j
	} else {
		h <- (j * .ABW_HASH_A + as.numeric(seed) * .ABW_HASH_B) %% .ABW_HASH_MOD
	}
	order(h)
}

#' Number of cross-fitting folds actually used
#'
#' @param n0 number of controls.
#' @param p number of covariate columns in \code{hbalobject$mat}.
#' @param nfolds requested number of folds (already validated).
#' @return \code{list(K = <integer>, message = <character or NULL>)}; the
#'   message is non-\code{NULL} only when \code{K} differs from \code{nfolds}.
#' @keywords internal
#' @noRd
.abw_choose_k <- function(n0, p, nfolds) {
	cap <- floor(n0 / (.ABW_CONTROLS_PER_FOLD_PER_COVARIATE * p))
	K <- as.integer(max(1, min(nfolds, cap)))
	msg <- NULL
	if (K != nfolds) {
		msg <- sprintf(
			"att(): using %d cross-fitting fold(s) (requested nfolds = %s) given %d controls and %d covariates.",
			K, format(nfolds, scientific = FALSE, trim = TRUE), as.integer(n0), as.integer(p))
	}
	list(K = K, message = msg)
}

#' Fold labels 1..K for the controls, in their row order
#'
#' @keywords internal
#' @noRd
.abw_make_folds <- function(n0, K, seed) {
	K <- as.integer(K)
	ord <- .abw_hash_order(n0, seed)
	fold <- integer(n0)
	fold[ord] <- ((seq_len(n0) - 1L) %% K) + 1L
	fold
}

#' Weighted ridge fit of the control outcome model on a standardized design
#'
#' Each column of \code{X} is centered and scaled with the training rows' own
#' (unweighted) mean and standard deviation; the intercept is never penalized.
#' With \code{lambda = 0} the solution is the weighted least-squares fit (it
#' reproduces \code{stats::lm.wfit()}'s fitted values on a full-rank design).
#'
#' @param X covariate matrix of the training rows (no intercept column).
#' @param Y outcomes of the training rows.
#' @param w weights of the training rows.
#' @param lambda ridge penalty (\code{>= 0}) on the standardized slopes.
#' @return \code{list(beta_std = <named numeric length ncol(X) + 1, in
#'   standardized units>, center = <numeric length ncol(X)>, scale = <numeric
#'   length ncol(X)>, rank = <integer, rank of the standardized design with
#'   intercept>)}.
#' @keywords internal
#' @noRd
.abw_ridge_fit <- function(X, Y, w, lambda) {
	p <- ncol(X)
	center <- colMeans(X)
	scl <- apply(X, 2, stats::sd)
	# Zero-variance (or, for a one-row design, undefined) columns are left
	# uncentered-scaled by 1; they then enter as all-zero columns. This should
	# not occur post-hbal (defensive only).
	scl[!is.finite(scl) | scl == 0] <- 1
	Xs <- scale(X, center = center, scale = scl)
	Zd <- cbind(1, Xs)
	# The rank is consumed only by the Hoerl-Kennard plug-in (full-control
	# sample, lambda = 0), as the residual degrees-of-freedom bookkeeping of
	# its variance estimate. It is not a validity gate: with lambda > 0 the
	# penalized normal equations are full rank whatever the rank of X.
	rank <- qr(Zd)$rank
	pen <- diag(c(0, rep(lambda, p)), nrow = p + 1L)
	# solve() can fail only for an exactly (or computationally) singular
	# system, i.e. lambda = 0 on a rank-deficient design. The fallback is the
	# intercept-only fit (predict the weighted training mean everywhere) --
	# numerical, not statistical, defensiveness.
	beta <- tryCatch(
		as.vector(solve(crossprod(Zd, Zd * w) + pen, crossprod(Zd, w * Y))),
		error = function(e) {
			sw <- sum(w)
			c(if (sw > 0) sum(w * Y) / sw else mean(Y), rep(0, p))
		}
	)
	cn <- colnames(X)
	if (is.null(cn)) cn <- paste0("X", seq_len(p))
	names(beta) <- c("(Intercept)", cn)
	list(beta_std = beta, center = center, scale = scl, rank = as.integer(rank))
}

#' Prediction from a ridge fit at new rows
#'
#' New rows are standardized with the fit's own training center/scale (never
#' their own mean/sd). There is no predict() method for this hand-rolled fit;
#' prediction is this direct matrix product.
#' @keywords internal
#' @noRd
.abw_ridge_predict <- function(fit, Xnew) {
	Xs <- scale(Xnew, center = fit$center, scale = fit$scale)
	as.vector(cbind(1, Xs) %*% fit$beta_std)
}

#' Ridge coefficients on the original scale of X
#'
#' Back-transforms the standardized coefficients so that
#' \code{cbind(1, X) \%*\% coef} equals \code{.abw_ridge_predict(fit, X)}.
#' @keywords internal
#' @noRd
.abw_coef_original <- function(fit) {
	b <- fit$beta_std
	slopes <- b[-1] / fit$scale
	out <- c(b[1] - sum(slopes * fit$center), slopes)
	names(out) <- names(b)
	out
}

#' Hoerl-Kennard ridge penalty from the full control sample
#'
#' Computed once per att() call, from the lambda = 0 (weighted least-squares)
#' fit on ALL controls: lambda = p * sigma2 / ||beta_ols[-1]||^2 with sigma2
#' the weighted residual variance, sum(w * resid^2) / max(sum(w) - rank, 1).
#' Falls back to lambda = 0 (plain weighted least squares) when sigma2 is not
#' finite and positive or when the slope norm is degenerate; capped at
#' \code{.ABW_LAMBDA_MAX}.
#'
#' @param X covariate matrix of all controls.
#' @param Y outcomes of all controls.
#' @param w weights of all controls.
#' @return \code{list(lambda = <numeric scalar >= 0>, rank = <integer>,
#'   sigma2 = <numeric>, beta_norm2 = <numeric>)}; only \code{lambda} is used
#'   by the estimator, the rest is diagnostic.
#' @keywords internal
#' @noRd
.abw_choose_lambda <- function(X, Y, w) {
	p <- ncol(X)
	fit0 <- .abw_ridge_fit(X, Y, w, lambda = 0)
	resid <- Y - .abw_ridge_predict(fit0, X)
	sigma2 <- sum(w * resid^2) / max(sum(w) - fit0$rank, 1)
	beta_norm2 <- sum(fit0$beta_std[-1]^2)
	lambda <- 0
	if (is.finite(sigma2) && sigma2 > 0 &&
	    is.finite(beta_norm2) && beta_norm2 > .ABW_BETA_NORM_TOL) {
		lambda <- max(0, min(p * sigma2 / beta_norm2, .ABW_LAMBDA_MAX))
	}
	list(lambda = lambda, rank = fit0$rank, sigma2 = sigma2, beta_norm2 = beta_norm2)
}

#' Cross-fitted augmented ATT estimator (the workhorse behind att(method = "abw"))
#'
#' @param hbalobject an \code{hbal} object that carries an outcome.
#' @param dr logical; \code{FALSE} drops the outcome model.
#' @param seed \code{NULL} or a single finite number (fold assignment only).
#' @param nfolds requested number of cross-fitting folds (already validated).
#' @return the ten-element list documented under \code{displayAll = TRUE} in
#'   \code{?att}: estimate, se, df, method, dr, influence, fold, nfolds,
#'   nuisance, weights.
#' @keywords internal
#' @noRd
.abw_att <- function(hbalobject, dr, seed, nfolds) {
	dr <- if (dr) TRUE else FALSE
	# The balancing weights enter for either dr value, so a non-converged
	# hbal() fit is flagged for both. hbal() stores $converged as an integer
	# 0/1, hence the as.logical() coercion (isTRUE(1L) is FALSE).
	if (!isTRUE(as.logical(hbalobject$converged))) {
		warning("att(): hbalobject$converged is FALSE -- the entropy-balancing weights did not converge. ",
		        "method = \"abw\"'s point estimate and standard error may be unreliable. Consider refitting ",
		        "hbal() with cv = TRUE, a lower expand.degree, or the ds / X.keep / exclude options.",
		        call. = FALSE)
	}
	Treatment <- hbalobject$Treatment
	tr <- Treatment == 1
	co <- Treatment == 0
	X <- as.matrix(hbalobject$mat)
	Y <- as.numeric(hbalobject$Outcome)
	w <- as.numeric(hbalobject$weights)
	gamma <- as.numeric(hbalobject$weights.co)
	n <- length(Treatment)
	n1 <- sum(tr)
	n0 <- sum(co)
	p <- ncol(X)
	if (nrow(X) != n || length(Y) != n || length(w) != n || length(gamma) != n0) {
		stop("att(): hbalobject$mat, $Outcome, $weights, and $weights.co have inconsistent lengths")
	}
	W1 <- sum(w[tr])
	X_tr <- X[tr, , drop = FALSE]
	X_co <- X[co, , drop = FALSE]
	Y_tr <- Y[tr]
	Y_co <- Y[co]
	w_tr <- w[tr]
	w_co <- w[co]

	fold <- rep(NA_integer_, n)
	if (dr) {
		kk <- .abw_choose_k(n0, p, nfolds)
		K <- kk$K
		if (!is.null(kk$message)) message(kk$message)
		fold_co <- .abw_make_folds(n0, K, seed)
		fold[co] <- fold_co
		# One ridge penalty for every fit, from the full control sample.
		lambda <- .abw_choose_lambda(X_co, Y_co, w_co)$lambda
		# m_i: outcome model fitted on all controls, evaluated at the treated rows
		fit_full <- .abw_ridge_fit(X_co, Y_co, w_co, lambda)
		m_hat <- .abw_ridge_predict(fit_full, X_tr)
		# e_i: cross-fitted control residuals (own-observation fit when K = 1)
		coef_folds <- vector("list", K)
		e_hat <- numeric(n0)
		if (K == 1L) {
			coef_folds[[1L]] <- .abw_coef_original(fit_full)
			e_hat <- Y_co - .abw_ridge_predict(fit_full, X_co)
		} else {
			for (k in seq_len(K)) {
				held <- fold_co == k
				fit_k <- .abw_ridge_fit(X_co[!held, , drop = FALSE], Y_co[!held], w_co[!held], lambda)
				coef_folds[[k]] <- .abw_coef_original(fit_k)
				e_hat[held] <- Y_co[held] - .abw_ridge_predict(fit_k, X_co[held, , drop = FALSE])
			}
		}
		nuisance <- list(coef_full = .abw_coef_original(fit_full), coef_folds = coef_folds,
		                 rank_full = fit_full$rank)
		K_used <- K
	} else {
		m_hat <- rep(0, n1)
		e_hat <- Y_co
		nuisance <- list(coef_full = NULL, coef_folds = list(), rank_full = NA_integer_)
		K_used <- NA_integer_
	}

	tau <- (sum(w_tr * (Y_tr - m_hat)) - sum(gamma * e_hat)) / W1
	psi <- numeric(n)
	psi[tr] <- w_tr * (Y_tr - m_hat - tau)
	psi[co] <- -gamma * e_hat
	V <- sum(psi^2) / W1^2

	list(
		estimate  = tau,
		se        = sqrt(V),
		df        = n - 1,
		method    = "abw",
		dr        = dr,
		influence = psi,
		fold      = fold,
		nfolds    = K_used,
		nuisance  = nuisance,
		weights   = list(treated = w_tr, control = w_co)
	)
}

#' The 1 x 7 display data.frame of att() built from the internal abw result
#'
#' @param res the list returned by \code{.abw_att()}.
#' @param Tr name of the treatment variable (row name of the result).
#' @keywords internal
#' @noRd
.abw_table <- function(res, Tr) {
	est <- res$estimate
	se <- res$se
	df <- res$df
	tval <- est / se
	crit <- stats::qt(.ABW_CI_PROB, df = df)
	out <- data.frame(
		"Estimate"   = est,
		"Std. Error" = se,
		"t value"    = tval,
		"Pr(>|t|)"   = 2 * stats::pt(-abs(tval), df = df),
		"CI Lower"   = est - crit * se,
		"CI Upper"   = est + crit * se,
		"DF"         = df,
		check.names = FALSE
	)
	rownames(out) <- Tr
	out
}
