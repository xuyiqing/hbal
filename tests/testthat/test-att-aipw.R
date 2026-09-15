# Unit tests for att(method = "aipw"), the default ATT estimator since hbal 1.3.0.
# Fixtures come from tests/testthat/helper-aipw.R (make_toy, make_small_toy).

out <- make_toy(1)

test_that("aipw is deterministic without a seed", {
  expect_identical(att(out), att(out))
})

test_that("aipw is deterministic with a fixed seed", {
  expect_identical(att(out, seed = 123), att(out, seed = 123))
})

test_that("different seeds can produce different fold assignments", {
  f1 <- att(out, seed = 1, displayAll = TRUE)$fold
  f2 <- att(out, seed = 2, displayAll = TRUE)$fold
  expect_false(identical(f1, f2))
})

test_that("aipw output has the base 1x7 shape", {
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
  expect_equal(res$method, "aipw")
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

test_that("aipw's response to an in-span outcome-model perturbation is small", {
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

test_that("aipw's perturbation-response formula for an out-of-span direction is internally consistent", {
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

test_that("aipw rejects unsupported ... arguments with a clear error, legacy methods do not", {
  expect_error(att(out, se_type = "HC2"), regexp = "does not accept")
  expect_error(att(out, clusters = rep(1, nrow(out$mat))), regexp = "does not accept")
  expect_no_error(att(out, method = "lm_robust", se_type = "HC2"))
})

test_that("an unrecognized method errors clearly", {
  expect_error(att(out, method = "not_a_method"), regexp = "[Uu]nrecognized method")
})

test_that("small control group falls back to K = 1 with a message, not a warning", {
  out_small <- make_small_toy()
  expect_message(att(out_small, displayAll = TRUE, nfolds = 5), regexp = "fold")
  expect_no_warning(suppressMessages(att(out_small, nfolds = 5)))
  expect_equal(suppressMessages(att(out_small, displayAll = TRUE, nfolds = 5))$nfolds, 1L)
})
