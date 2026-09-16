# Unit tests for att(method = "abw"), the default ATT estimator since hbal 1.3.0.
# Fixtures come from tests/testthat/helper-abw.R (make_toy, make_small_toy,
# make_wide_toy).

out <- make_toy(1)

test_that("abw is deterministic without a seed", {
  expect_identical(att(out), att(out))
})

test_that("abw is deterministic with a fixed seed", {
  expect_identical(att(out, seed = 123), att(out, seed = 123))
})

test_that("different seeds can produce different fold assignments", {
  f1 <- att(out, seed = 1, displayAll = TRUE)$fold
  f2 <- att(out, seed = 2, displayAll = TRUE)$fold
  expect_false(identical(f1, f2))
})

test_that("abw output has the base 1x7 shape", {
  res <- att(out)
  expect_s3_class(res, "data.frame")
  expect_false(inherits(res, "tbl"))
  expect_equal(dim(res), c(1, 7))
  expect_named(res, c("Estimate", "Std. Error", "t value", "Pr(>|t|)",
                      "CI Lower", "CI Upper", "DF"))
  expect_equal(rownames(res), out$Treat)
  expect_true(all(is.finite(unlist(res))))
})

test_that("displayAll = TRUE returns the documented list shape", {
  res <- att(out, displayAll = TRUE)
  expect_type(res, "list")
  expect_false(is.data.frame(res))
  expect_named(res, c("estimate", "se", "df", "method", "dr", "influence",
                      "fold", "nfolds", "nuisance", "weights"),
               ignore.order = FALSE)
  expect_length(res$influence, nrow(out$mat))
  expect_named(res$nuisance, c("coef_full", "coef_folds", "rank_full"))
  expect_length(res$weights$treated, sum(out$Treatment == 1))
  expect_length(res$weights$control, sum(out$Treatment == 0))
  expect_equal(res$method, "abw")
})

test_that("dr = FALSE reduces to the weighted difference in means", {
  ref <- (sum(out$weights[out$Treatment == 1] * out$Outcome[out$Treatment == 1]) -
          sum(out$weights[out$Treatment == 0] * out$Outcome[out$Treatment == 0])) /
         sum(out$weights[out$Treatment == 1])
  expect_equal(att(out, dr = FALSE)$Estimate, ref, tolerance = 1e-8)
})

test_that("closed-form agreement under exact balance and no cross-fitting", {
  # n0 = 10 < 4p = 12 guarantees the K = 1 fallback, so the exact-balance
  # collapse of the augmented estimator onto the weighted difference in means
  # applies without cross-fitting slack. Tolerance 0.02 accommodates hbal's own
  # solver slack around constraint.tolerance = 1e-3, not model error.
  out_small <- make_small_toy()
  # hbal() stores $converged as the integer 1/0 returned by the C++ solver, so
  # the guard coerces before isTRUE() (a bare isTRUE(1L) is FALSE).
  testthat::skip_if_not(isTRUE(as.logical(out_small$converged)),
                        "hbal did not converge; closed-form check does not apply")
  ref <- (sum(out_small$weights[out_small$Treatment == 1] * out_small$Outcome[out_small$Treatment == 1]) -
          sum(out_small$weights[out_small$Treatment == 0] * out_small$Outcome[out_small$Treatment == 0])) /
         sum(out_small$weights[out_small$Treatment == 1])
  expect_equal(att(out_small, dr = TRUE)$Estimate, ref, tolerance = 0.02)
})

test_that("abw's response to an in-span outcome-model perturbation is small", {
  # White-box check of the new estimator only: perturb the fitted nuisance
  # PREDICTIONS by eps * h (nothing is refit, $Outcome is untouched). For a
  # direction h inside the exactly balanced linear block, the moment formula
  # gives delta_tau = -eps * (weighted treated mean of h - weighted control
  # mean of h), which hbal's balance makes small.
  out <- make_toy(2026, n1 = 300, n0 = 300, expand.degree = 1, cv = FALSE)
  res <- att(out, displayAll = TRUE)
  h <- out$mat[, 1]   # first column of $mat: always inside the exactly balanced linear block
  eps <- 1
  imbalance <- sum(res$weights$treated * h[out$Treatment == 1]) / sum(res$weights$treated) -
               sum(res$weights$control * h[out$Treatment == 0]) / sum(res$weights$control)
  delta_tau <- -eps * imbalance
  expect_lt(abs(delta_tau), 0.02)
})

test_that("abw's perturbation-response formula for an out-of-span direction is internally consistent", {
  # Fully independent recomputation (nothing reused from the previous test):
  # displayAll = TRUE's exposed weights must agree exactly with hbalobject$weights,
  # so the two representations of "the same weight" cannot silently diverge.
  # This does not claim orthogonality for an out-of-span direction.
  out <- make_toy(2026, n1 = 300, n0 = 300, expand.degree = 1, cv = FALSE)
  res <- att(out, displayAll = TRUE)
  h <- out$mat[, 1] * out$mat[, 2]   # NOT a column of out$mat when expand.degree = 1 (deliberately)
  eps <- 1
  imbalance_from_displayAll <- sum(res$weights$treated * h[out$Treatment == 1]) / sum(res$weights$treated) -
                               sum(res$weights$control * h[out$Treatment == 0]) / sum(res$weights$control)
  w_tr <- out$weights[out$Treatment == 1]; w_co <- out$weights[out$Treatment == 0]
  h_tr <- h[out$Treatment == 1];           h_co <- h[out$Treatment == 0]
  imbalance_from_raw_object <- sum(w_tr * h_tr) / sum(w_tr) - sum(w_co * h_co) / sum(w_co)
  expect_equal(imbalance_from_displayAll, imbalance_from_raw_object, tolerance = 1e-10)
  expect_true(is.finite(-eps * imbalance_from_raw_object))
})

test_that("legacy lm_robust point estimate matches an independently-computed weighted OLS fit", {
  # Valid because lm_robust's se_type only changes the variance formula, never
  # the WLS point estimate (same normal equations regardless of se_type).
  dat <- cbind.data.frame(out$Outcome, out$Treatment, out$mat)
  colnames(dat) <- c(out$Y, out$Treat, colnames(out$mat))
  ff <- as.formula(paste0(out$Y, " ~ ", out$Treat, " + ", paste(colnames(out$mat), collapse = " + ")))
  ref <- unname(coef(lm(ff, data = dat, weights = out$weights))[out$Treat])
  expect_equal(att(out, method = "lm_robust", dr = TRUE)$Estimate, ref, tolerance = 1e-8)
})

test_that("legacy lm_lin runs and returns the base shape", {
  res <- att(out, method = "lm_lin")
  expect_s3_class(res, "data.frame")
  expect_false(inherits(res, "tbl"))
  expect_equal(dim(res), c(1, 7))
  expect_named(res, c("Estimate", "Std. Error", "t value", "Pr(>|t|)",
                      "CI Lower", "CI Upper", "DF"))
  expect_equal(rownames(res), out$Treat)
  expect_true(all(is.finite(unlist(res))))
})

test_that("elnet is reproducible when the caller fixes the seed", {
  # set.seed() here is caller/test code, not package code.
  set.seed(999); r1 <- att(out, method = "elnet")
  set.seed(999); r2 <- att(out, method = "elnet")
  expect_equal(r1, r2)
})

test_that("abw rejects unsupported ... arguments with a clear error, legacy methods do not", {
  expect_error(att(out, se_type = "HC2"), regexp = "does not accept")
  expect_error(att(out, clusters = rep(1, nrow(out$mat))), regexp = "does not accept")
  expect_no_error(att(out, method = "lm_robust", se_type = "HC2"))
})

test_that("an unrecognized method errors clearly", {
  expect_error(att(out, method = "not_a_method"), regexp = "[Uu]nrecognized method")
})

test_that("small control group falls back to K = 1 with a message, not a warning", {
  out_small <- make_small_toy()
  # The "no warning" claim below is only meaningful on a converged fit (a
  # non-converged one warns by design, see the next test).
  expect_true(isTRUE(as.logical(out_small$converged)))
  expect_message(att(out_small, displayAll = TRUE, nfolds = 5), regexp = "fold")
  expect_no_warning(suppressMessages(att(out_small, nfolds = 5)))
  expect_equal(suppressMessages(att(out_small, displayAll = TRUE, nfolds = 5))$nfolds, 1L)
})

test_that("a non-converged hbalobject triggers exactly the documented warning, for both dr values", {
  out_nonconv <- make_wide_toy()
  testthat::skip_if(isTRUE(as.logical(out_nonconv$converged)),
                    "the wide fixture converged on this platform; the non-convergence path cannot be exercised")
  # suppressMessages(): n0 = 60, p = 15 also triggers the (separate) fold-count message.
  expect_warning(suppressMessages(att(out_nonconv)), regexp = "converged")
  expect_warning(att(out_nonconv, dr = FALSE), regexp = "converged")
  expect_length(testthat::capture_warnings(suppressMessages(att(out_nonconv))), 1L)
  expect_length(testthat::capture_warnings(att(out_nonconv, dr = FALSE)), 1L)
  # and a converged fit stays silent, for both dr values
  expect_no_warning(att(out))
  expect_no_warning(att(out, dr = FALSE))
})

test_that("abw's point estimate stays bounded on a wide/small-effective-sample control design", {
  # Coarse regression guard tied to the outcome's own scale, not a tight
  # statistical claim. On this fixture a plain weighted least-squares outcome
  # model (ridge penalty 0) gives an ATT of about 44 x range(Y); the ridge
  # outcome model with the GCV-selected penalty gives about 0.3 x range(Y).
  out_nonconv <- make_wide_toy()
  testthat::skip_if(isTRUE(as.logical(out_nonconv$converged)),
                    "the wide fixture converged on this platform; it does not stress the intended mechanism")
  est <- suppressWarnings(suppressMessages(att(out_nonconv)))$Estimate
  expect_true(is.finite(est))
  expect_lt(abs(est), 20 * diff(range(out_nonconv$Outcome)))
})

test_that("the ridge outcome model at lambda = 0 reproduces lm.wfit() on a full-rank design", {
  co <- out$Treatment == 0
  X <- as.matrix(out$mat)[co, , drop = FALSE]
  Y <- out$Outcome[co]
  w <- out$weights[co]
  fit <- hbal:::.abw_ridge_fit(X, Y, w, lambda = 0)
  ref <- stats::lm.wfit(x = cbind(1, X), y = Y, w = w)
  expect_equal(fit$rank, ref$rank)
  expect_equal(unname(hbal:::.abw_ridge_predict(fit, X)), unname(ref$fitted.values), tolerance = 1e-10)
  # the original-scale coefficients reproduce the same predictions
  expect_equal(unname(as.vector(cbind(1, X) %*% hbal:::.abw_coef_original(fit))),
               unname(ref$fitted.values), tolerance = 1e-10)
  # and a positive penalty shrinks the standardized slopes but not the intercept's role
  fit5 <- hbal:::.abw_ridge_fit(X, Y, w, lambda = 5)
  expect_lt(sum(fit5$beta_std[-1]^2), sum(fit$beta_std[-1]^2))
  expect_true(all(is.finite(hbal:::.abw_ridge_predict(fit5, X))))
})

test_that("nuisance$coef_full is on the scale of mat and reproduces the treated-unit predictions implied by influence", {
  res <- att(out, displayAll = TRUE)
  tr <- out$Treatment == 1
  # psi_i = w_i (Y_i - m_i - tau) for treated rows  =>  m_i = Y_i - psi_i / w_i - tau
  m_from_influence <- out$Outcome[tr] - res$influence[tr] / res$weights$treated - res$estimate
  m_from_coef <- as.vector(cbind(1, as.matrix(out$mat)[tr, , drop = FALSE]) %*% res$nuisance$coef_full)
  expect_equal(m_from_coef, m_from_influence, tolerance = 1e-10)
  expect_named(res$nuisance$coef_full, c("(Intercept)", colnames(out$mat)))
  expect_length(res$nuisance$coef_folds, res$nfolds)
})

test_that("abw's ridge outcome model stays bounded at the exact interpolation boundary (regression: tester's BLOCK-2 fixture)", {
  # Verbatim reproduction of the fixture on which a plug-in ridge penalty
  # collapsed to ~0: the full-control design sits at the exact interpolation
  # boundary (rank_full == n0 == 15), so the unpenalized fit interpolates the
  # controls and extrapolates wildly to the treated rows (att(out)$Estimate was
  # -3409.35, about 320 x diff(range(Outcome))). With the penalty chosen by GCV
  # under the effective-degrees-of-freedom cap, such a fit is inadmissible.
  set.seed(5)
  n1 <- 80; n0 <- 15; p_raw <- 4
  n <- n1 + n0
  Xs <- as.data.frame(matrix(rnorm(n * p_raw), n, p_raw))
  names(Xs) <- paste0("X", 1:p_raw)
  Tr <- c(rep(1, n1), rep(0, n0))
  Y <- 1 + 0.5 * Tr + rowSums(Xs) + rnorm(n)
  dat <- data.frame(Tr = Tr, Xs, Y = Y)
  out <- suppressMessages(hbal::hbal(Treat = "Tr", X = names(Xs), Y = "Y", data = dat,
                                     expand.degree = 2, print.level = -1))
  expect_equal(ncol(out$mat), 14L)          # documents the fixture: rank_full == n0 == 15 exactly
  # suppressWarnings(): a non-convergence warning is expected and irrelevant here;
  # suppressMessages(): n0 = 15 < 2p also triggers the K = 1 fold-count message.
  res <- suppressWarnings(suppressMessages(att(out)))
  expect_true(all(is.finite(unlist(res))))
  expect_lt(abs(res$Estimate), 20 * diff(range(out$Outcome)))   # the bounded-estimate invariant
  expect_lt(abs(res$Estimate), 1 * diff(range(out$Outcome)))    # much tighter, evidence-based: the
                                                                 # GCV rule gives a ratio of about 0.075
                                                                 # on this exact fixture; this catches a
                                                                 # regression toward "technically under
                                                                 # 20x but still huge" long before the
                                                                 # 20x floor would
  full <- suppressWarnings(suppressMessages(att(out, displayAll = TRUE)))
  expect_equal(full$nuisance$rank_full, 15L)   # the boundary itself: rank_full == n0
})

test_that("the new lambda rule leaves well-conditioned estimates close to the previous plug-in's", {
  # n0 = 150 >> 2p = 6: the effective-degrees-of-freedom cap never binds here,
  # and GCV should agree closely with the value the Hoerl-Kennard plug-in of
  # the previous revision gave on this exact fixture (0.5082179957). Guards
  # against either rule over-shrinking ordinary, well-conditioned data.
  res <- att(out)
  expect_equal(res$Estimate, 0.5082179957, tolerance = 0.01)
  lam <- hbal:::.abw_choose_lambda(as.matrix(out$mat)[out$Treatment == 0, , drop = FALSE],
                                    out$Outcome[out$Treatment == 0],
                                    out$weights[out$Treatment == 0])
  expect_lt(lam$lambda, 10)
  # the diagnostic helper's documented return contract
  expect_named(lam, c("lambda", "edf", "cap", "n_eligible"))
  expect_true(lam$lambda %in% hbal:::.ABW_LAMBDA_GRID)
  expect_lte(lam$edf, lam$cap + 1e-8)
  expect_equal(lam$cap, min(ncol(out$mat), sum(out$Treatment == 0) / 2) + 1)
})
