# Internal helpers for att(method = "abw"): the cross-fitted, Neyman-orthogonal
# augmented balancing weights (ABW) ATT estimator that is the default of att()
# since hbal 1.3.0. Nothing in this file is exported, and no function here
# reads or writes .Random.seed (no set.seed(), sample(), runif(), ...): the fold
# assignment is a deterministic hash of the `seed` argument, and the ridge
# penalty of the outcome model is chosen by a deterministic grid search (GCV;
# no data splitting).
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
# lambda >= 0 is computed ONCE per att() call from the full control sample and
# reused for every fit: generalized cross-validation (GCV; Golub, Heath and
# Wahba 1979) minimized over the fixed grid .ABW_LAMBDA_GRID, restricted to
# the candidates whose effective degrees of freedom
#   edf(lambda) = tr[(A + pen)^{-1} A],   A = Zd' W Zd
# do not exceed min(p, n0 / 2) + 1 (the + 1 is the unpenalized intercept), so
# the fit can never approach an interpolation of the controls. No data
# splitting and no random numbers are involved. lambda = 0 reproduces weighted
# least squares exactly.

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
# Candidate ridge penalties (standardized design) searched by GCV: 0 and
# 10^j for j = -8, -7.75, ..., 6 (58 points in all). The grid's upper bound,
# 1e6, is the most regularization the outcome model can receive.
.ABW_LAMBDA_GRID <- c(0, 10^seq(-8, 6, by = 0.25))
# Floating-point slack on the effective-degrees-of-freedom cap: edf(0) equals
# p + 1 in exact arithmetic but is computed by a linear solve, and the cap is
# exactly p + 1 whenever n0 >= 2p, where the unpenalized fit must stay
# admissible; the comparison edf <= cap therefore gets this much tolerance.
.ABW_EDF_TOL <- 1e-8
# The GCV score n0 * RSS / (n0 - edf)^2 counts as +Inf when its denominator
# is not larger than this.
.ABW_GCV_DENOM_TOL <- 1e-8

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

#' Standardize the columns of a design matrix with their own mean and sd
#'
#' The one standardization used in this file, by the ridge fit (on its
#' training rows) and by the penalty choice (on the full control sample), so
#' that the design GCV evaluates is the design the fit uses.
#'
#' @param X covariate matrix (no intercept column).
#' @return \code{list(Xs = <standardized matrix>, center = <numeric length
#'   ncol(X)>, scale = <numeric length ncol(X)>)}.
#' @keywords internal
#' @noRd
.abw_standardize <- function(X) {
	center <- colMeans(X)
	scl <- apply(X, 2, stats::sd)
	# Zero-variance (or, for a one-row design, undefined) columns get scale 1;
	# once centered they enter as all-zero columns. This should not occur
	# post-hbal (defensive only).
	scl[!is.finite(scl) | scl == 0] <- 1
	list(Xs = scale(X, center = center, scale = scl), center = center, scale = scl)
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
	std <- .abw_standardize(X)
	center <- std$center
	scl <- std$scale
	Zd <- cbind(1, std$Xs)
	# The rank of the standardized design (with intercept) is a diagnostic
	# only, exposed as nuisance$rank_full for the full-control fit. It is not
	# a validity gate and plays no part in the penalty choice: with lambda > 0
	# the penalized normal equations are full rank whatever the rank of X.
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

#' Ridge penalty by GCV on a fixed grid, gated by an effective-degrees-of-freedom cap
#'
#' Computed once per att() call from ALL controls, on the same standardized
#' design the ridge fit uses (\code{.abw_standardize()}). For every candidate
#' in \code{.ABW_LAMBDA_GRID} the effective degrees of freedom of the ridge
#' hat matrix, \code{edf(lambda) = tr[(A + pen)^-1 A]} with
#' \code{A = Zd' W Zd}, is computed; a candidate is admissible when
#' \code{edf(lambda) <= min(p, n0 / 2) + 1} (the + 1 is the never-penalized
#' intercept, so the unpenalized fit, edf = p + 1, stays admissible whenever
#' n0 >= 2p). Among the admissible candidates the Golub, Heath and Wahba
#' (1979) GCV score \code{n0 * RSS(lambda) / (n0 - edf(lambda))^2} is
#' minimized; ties go to the smallest lambda (the grid is ascending and
#' \code{which.min()} returns the first minimum). Fallbacks, none of which
#' can occur on a valid hbal object (hbal() guarantees n0 >= p + 1): a
#' candidate whose linear solve fails or whose edf is not finite is
#' inadmissible; an empty admissible set falls back to the largest grid
#' value; a candidate whose GCV denominator is not positive, or whose score
#' is not finite, scores +Inf; when every admissible score is +Inf the
#' largest admissible lambda is used.
#'
#' @param X covariate matrix of all controls.
#' @param Y outcomes of all controls.
#' @param w weights of all controls.
#' @return \code{list(lambda = <numeric scalar >= 0>, edf = <numeric, the
#'   effective degrees of freedom at lambda>, cap = <numeric, the edf cap>,
#'   n_eligible = <integer, number of admissible grid points>)}; only
#'   \code{lambda} is used by the estimator, the rest is diagnostic.
#' @keywords internal
#' @noRd
.abw_choose_lambda <- function(X, Y, w) {
	p <- ncol(X)
	n0 <- nrow(X)
	Zd <- cbind(1, .abw_standardize(X)$Xs)
	A <- crossprod(Zd, Zd * w)
	b <- crossprod(Zd, w * Y)
	penalty <- function(lambda) diag(c(0, rep(lambda, p)), nrow = p + 1L)
	# Effective degrees of freedom of the ridge hat matrix: p + 1 at lambda = 0
	# (unpenalized) down to 1 as lambda -> Inf (intercept only). A failed solve
	# (lambda = 0 on a rank-deficient design) makes the candidate inadmissible
	# rather than an error.
	edf <- function(lambda) {
		tryCatch(sum(diag(solve(A + penalty(lambda), A))), error = function(e) NA_real_)
	}
	cap <- min(p, n0 / 2) + 1
	edf_grid <- vapply(.ABW_LAMBDA_GRID, edf, numeric(1))
	admissible <- is.finite(edf_grid) & edf_grid <= cap + .ABW_EDF_TOL
	eligible <- .ABW_LAMBDA_GRID[admissible]
	edf_eligible <- edf_grid[admissible]
	if (length(eligible) == 0L) {
		eligible <- max(.ABW_LAMBDA_GRID)
		edf_eligible <- edf(eligible)
	}
	gcv_score <- function(lambda, edf_lambda) {
		beta <- tryCatch(solve(A + penalty(lambda), b), error = function(e) NULL)
		if (is.null(beta)) return(Inf)
		rss <- sum(w * (Y - Zd %*% beta)^2)
		denom <- n0 - edf_lambda
		if (is.finite(denom) && denom > .ABW_GCV_DENOM_TOL) n0 * rss / denom^2 else Inf
	}
	scores <- vapply(seq_along(eligible),
	                 function(i) gcv_score(eligible[i], edf_eligible[i]), numeric(1))
	scores[!is.finite(scores)] <- Inf
	i <- if (all(is.infinite(scores))) length(eligible) else which.min(scores)
	list(lambda = eligible[i], edf = edf_eligible[i], cap = cap,
	     n_eligible = length(eligible))
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
