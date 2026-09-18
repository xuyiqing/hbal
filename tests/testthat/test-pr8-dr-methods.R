# Regression test for GitHub PR #8 (soodoku): att(dr = FALSE) was accepted and
# silently ignored by method = "lm_lin" and method = "elnet", and the guard that
# was meant to catch it compared method against the unused string "lin".

test_that("dr = FALSE gives the weighted difference in means for every method", {
	set.seed(3)
	n <- 400
	x1 <- rnorm(n)
	x2 <- rnorm(n)
	treat <- rbinom(n, 1, plogis(0.8 * x1))
	Y <- 0.5 * treat + x1 + 0.5 * x2 + rnorm(n)
	dat <- data.frame(treat = treat, Y = Y, cov1 = x1, cov2 = x2)
	out <- hbal(Treat = "treat", X = c("cov1", "cov2"), Y = "Y", data = dat,
				print.level = -1)

	tr <- out$Treatment == 1
	w <- out$weights
	ref <- weighted.mean(out$Outcome[tr], w[tr]) -
		weighted.mean(out$Outcome[!tr], w[!tr])

	base <- att(out, method = "lm_robust", dr = FALSE)
	expect_equal(base$Estimate, ref, tolerance = 1e-8)

	for (m in c("lm_lin", "elnet")) {
		res <- att(out, method = m, dr = FALSE)
		expect_equal(res$Estimate, ref, tolerance = 1e-8)
		expect_equal(as.numeric(res), as.numeric(base), tolerance = 1e-10)
	}
})

test_that("dr = FALSE leaves the abw default on its own path", {
	set.seed(3)
	n <- 400
	x1 <- rnorm(n)
	x2 <- rnorm(n)
	treat <- rbinom(n, 1, plogis(0.8 * x1))
	Y <- 0.5 * treat + x1 + 0.5 * x2 + rnorm(n)
	dat <- data.frame(treat = treat, Y = Y, cov1 = x1, cov2 = x2)
	out <- hbal(Treat = "treat", X = c("cov1", "cov2"), Y = "Y", data = dat,
				print.level = -1)

	base <- att(out, method = "lm_robust", dr = FALSE)
	abw <- att(out, dr = FALSE)
	expect_equal(abw$Estimate, base$Estimate, tolerance = 1e-8)
	expect_false(isTRUE(all.equal(abw$DF, base$DF)))
	expect_false(isTRUE(all.equal(abw[["Std. Error"]], base[["Std. Error"]])))

	full <- att(out, dr = FALSE, displayAll = TRUE)
	expect_false(inherits(full, "lm_robust"))
	expect_identical(full$method, "abw")
})

test_that("dr = TRUE is unchanged for lm_lin", {
	set.seed(3)
	n <- 400
	x1 <- rnorm(n)
	x2 <- rnorm(n)
	treat <- rbinom(n, 1, plogis(0.8 * x1))
	Y <- 0.5 * treat + x1 + 0.5 * x2 + rnorm(n)
	dat <- data.frame(treat = treat, Y = Y, cov1 = x1, cov2 = x2)
	out <- hbal(Treat = "treat", X = c("cov1", "cov2"), Y = "Y", data = dat,
				print.level = -1)

	fit <- att(out, method = "lm_lin", dr = TRUE, displayAll = TRUE)
	expect_true(any(grepl("^treat:", names(stats::coef(fit)))))
	expect_false(isTRUE(all.equal(
		as.numeric(att(out, method = "lm_lin", dr = TRUE)),
		as.numeric(att(out, method = "lm_lin", dr = FALSE)))))
})

test_that('method = "lin" is rejected whatever dr is', {
	set.seed(3)
	n <- 200
	x1 <- rnorm(n)
	treat <- rbinom(n, 1, plogis(0.8 * x1))
	Y <- 0.5 * treat + x1 + rnorm(n)
	dat <- data.frame(treat = treat, Y = Y, cov1 = x1)
	out <- hbal(Treat = "treat", X = "cov1", Y = "Y", data = dat, print.level = -1)

	expect_error(att(out, method = "lin", dr = TRUE), "unrecognized method")
	expect_error(att(out, method = "lin", dr = FALSE), "unrecognized method")
})
