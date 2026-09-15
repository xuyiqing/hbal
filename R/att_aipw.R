# Internal helpers for att(method = "aipw"): the cross-fitted, Neyman-orthogonal
# augmented balancing-weights (AIPW-style) ATT estimator that is the default of
# att() since hbal 1.3.0. Nothing in this file is exported, and no function here
# reads or writes .Random.seed (no set.seed(), sample(), runif(), ...): the fold
# assignment is a deterministic hash of the `seed` argument.
#
# Notation (i indexes the rows of the hbal object in their original order):
#   T_i      treatment indicator                        hbalobject$Treatment
#   X_i      covariate row, unscaled, already expanded   hbalobject$mat
#   Y_i      outcome                                     hbalobject$Outcome
#   w_i      base weight of every unit                   hbalobject$weights
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

# Multiplicative-hash constants (Knuth / xxHash primes) used to turn `seed`
# into a deterministic permutation of the controls without touching R's RNG.
.AIPW_HASH_A <- 2654435761
.AIPW_HASH_B <- 2246822519
.AIPW_HASH_MOD <- 4294967296  # 2^32
# Upper quantile of the t distribution for the two-sided 95 percent interval.
.AIPW_CI_PROB <- 0.975
# Cross-fitting keeps at least this many controls per fold per covariate
# column, so that fold-level fits stay reasonably determined:
# K <= floor(n0 / (.AIPW_CONTROLS_PER_FOLD_PER_COVARIATE * p)).
.AIPW_CONTROLS_PER_FOLD_PER_COVARIATE <- 2

#' Deterministic permutation of the controls used for fold assignment
#'
#' @param n0 number of controls.
#' @param seed \code{NULL} (identity order) or a single finite number.
#' @return an integer permutation of \code{seq_len(n0)}.
#' @keywords internal
#' @noRd
.aipw_hash_order <- function(n0, seed) {
	j <- as.numeric(seq_len(n0))
	if (is.null(seed)) {
		h <- j
	} else {
		h <- (j * .AIPW_HASH_A + as.numeric(seed) * .AIPW_HASH_B) %% .AIPW_HASH_MOD
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
.aipw_choose_k <- function(n0, p, nfolds) {
	cap <- floor(n0 / (.AIPW_CONTROLS_PER_FOLD_PER_COVARIATE * p))
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
.aipw_make_folds <- function(n0, K, seed) {
	K <- as.integer(K)
	ord <- .aipw_hash_order(n0, seed)
	fold <- integer(n0)
	fold[ord] <- ((seq_len(n0) - 1L) %% K) + 1L
	fold
}

#' Weighted least-squares fit of the control outcome model
#'
#' @param X covariate matrix of the training rows (no intercept column).
#' @param Y outcomes of the training rows.
#' @param w base weights of the training rows.
#' @return \code{list(coef = <numeric length ncol(X) + 1, NA-free>, rank = <integer>)}.
#' @keywords internal
#' @noRd
.aipw_fit <- function(X, Y, w) {
	# Rank deficiency (aliased columns of a wide or collinear expanded design,
	# or a fold-level training set too small to identify every column) is
	# handled entirely by lm.wfit()'s pivoted QR: aliased coefficients come
	# back as NA, already un-pivoted into the ORIGINAL column order, and are
	# zeroed here, so that cbind(1, X) %*% coef is the least-squares prediction
	# from the retained columns. No explicit rank check is made and no `tol` is
	# passed (R's default 1e-7 applies). .lm.fit() must NOT be substituted
	# here: its coefficients come back in pivoted order.
	fit <- stats::lm.wfit(x = cbind(1, X), y = Y, w = w)
	beta <- fit$coefficients
	beta[is.na(beta)] <- 0
	cn <- colnames(X)
	if (is.null(cn)) cn <- paste0("X", seq_len(ncol(X)))
	names(beta) <- c("(Intercept)", cn)
	list(coef = beta, rank = as.integer(fit$rank))
}

#' Linear prediction cbind(1, X) times coef
#'
#' lm.wfit() objects have no predict() method, so prediction is this direct
#' matrix product.
#' @keywords internal
#' @noRd
.aipw_predict <- function(coef, X) {
	as.vector(cbind(1, X) %*% coef)
}

#' Cross-fitted augmented ATT estimator (the workhorse behind att(method = "aipw"))
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
.aipw_att <- function(hbalobject, dr, seed, nfolds) {
	dr <- if (dr) TRUE else FALSE
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
		kk <- .aipw_choose_k(n0, p, nfolds)
		K <- kk$K
		if (!is.null(kk$message)) message(kk$message)
		fold_co <- .aipw_make_folds(n0, K, seed)
		fold[co] <- fold_co
		# m_i: outcome model fitted on all controls, evaluated at the treated rows
		fit_full <- .aipw_fit(X_co, Y_co, w_co)
		m_hat <- .aipw_predict(fit_full$coef, X_tr)
		# e_i: cross-fitted control residuals (own-observation fit when K = 1)
		coef_folds <- vector("list", K)
		e_hat <- numeric(n0)
		if (K == 1L) {
			coef_folds[[1L]] <- fit_full$coef
			e_hat <- Y_co - .aipw_predict(fit_full$coef, X_co)
		} else {
			for (k in seq_len(K)) {
				held <- fold_co == k
				fit_k <- .aipw_fit(X_co[!held, , drop = FALSE], Y_co[!held], w_co[!held])
				coef_folds[[k]] <- fit_k$coef
				e_hat[held] <- Y_co[held] - .aipw_predict(fit_k$coef, X_co[held, , drop = FALSE])
			}
		}
		nuisance <- list(coef_full = fit_full$coef, coef_folds = coef_folds, rank_full = fit_full$rank)
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
		method    = "aipw",
		dr        = dr,
		influence = psi,
		fold      = fold,
		nfolds    = K_used,
		nuisance  = nuisance,
		weights   = list(treated = w_tr, control = w_co)
	)
}

#' The 1 x 7 display data.frame of att() built from the internal aipw result
#'
#' @param res the list returned by \code{.aipw_att()}.
#' @param Tr name of the treatment variable (row name of the result).
#' @keywords internal
#' @noRd
.aipw_table <- function(res, Tr) {
	est <- res$estimate
	se <- res$se
	df <- res$df
	tval <- est / se
	crit <- stats::qt(.AIPW_CI_PROB, df = df)
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
