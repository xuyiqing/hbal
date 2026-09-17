# Regression test for GitHub PR #7 (soodoku). hbal() renames the covariates in X
# to X1..Xk internally. Two defects in that one statement made it balance
# something other than what the caller asked for: the generated names could
# duplicate a column that was NOT being renamed, and the names were assigned in
# data-frame order while everything downstream assumes user X order.

dev_co <- function(fit, d, v) {
	w <- fit$weights.co
	co <- d$treat == 0
	tr <- d$treat == 1
	abs(mean(d[[v]][tr]) - sum(w * d[[v]][co]) / sum(w))
}

test_that("a pre-existing X1 column does not displace the requested covariates", {
	set.seed(11)
	n <- 600
	X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n); X4 <- rnorm(n)
	treat <- rbinom(n, 1, plogis(0.9 * X2 + 0.7 * X3))
	Y <- 0.5 * treat + X2 + X3 + rnorm(n)
	dat <- data.frame(X1, X2, X3, X4, treat, Y)

	full <- hbal(Treat = "treat", X = c("X2", "X3"), Y = "Y", data = dat,
				 print.level = -1)
	sub <- hbal(Treat = "treat", X = c("X2", "X3"), Y = "Y",
				data = dat[, c("treat", "X2", "X3", "Y")], print.level = -1)

	expect_lt(dev_co(full, dat, "X2"), 1e-3)
	expect_lt(dev_co(full, dat, "X3"), 1e-3)
	expect_equal(unname(full$weights.co), unname(sub$weights.co), tolerance = 1e-6)
	expect_equal(att(full, method = "lm_robust")$Estimate,
				 att(sub, method = "lm_robust")$Estimate, tolerance = 1e-6)
	expect_identical(rownames(full$bal.tab), c("X2", "X3"))
})

test_that("covariates are labelled in the order the caller supplied them", {
	set.seed(21)
	n <- 400
	a <- rnorm(n); b <- rnorm(n); cc <- rnorm(n)
	treat <- rbinom(n, 1, plogis(0.6 * a + 0.5 * cc))
	Y <- 0.5 * treat + a + cc + rnorm(n)
	dat <- data.frame(a, b, cc, treat, Y)

	fit <- hbal(Treat = "treat", X = c("cc", "a", "b"), Y = "Y", data = dat,
				print.level = -1)
	expect_identical(rownames(fit$bal.tab), c("cc", "a", "b"))
	expect_equal(unname(fit$bal.tab[, "Tr.Mean"]),
				 round(c(mean(cc[treat == 1]), mean(a[treat == 1]),
						 mean(b[treat == 1])), 2))
})

test_that("X.expand expands the covariate that was named", {
	set.seed(31)
	n <- 300
	a <- rnorm(n); b <- rnorm(n)
	treat <- rbinom(n, 1, plogis(0.5 * a))
	Y <- 0.5 * treat + a + b + rnorm(n)
	dat <- data.frame(a, b, treat, Y)

	fit <- hbal(Treat = "treat", X = c("b", "a"), Y = "Y", data = dat,
				expand.degree = 2, X.expand = "b", print.level = -1)
	expect_equal(unname(fit$mat[, 3]), b^2, tolerance = 1e-8)
})

test_that("the generated prefix never reaches user-facing names", {
	set.seed(61)
	n <- 500
	X1 <- rnorm(n); p <- rnorm(n); q <- rnorm(n)
	treat <- rbinom(n, 1, plogis(0.7 * p))
	Y <- 0.5 * treat + p + q + rnorm(n)
	dat <- data.frame(X1 = X1, p = p, q = q, treat = treat, Y = Y)

	fit <- hbal(Treat = "treat", X = c("p", "q"), Y = "Y", data = dat,
				expand.degree = 2, print.level = -1)
	expect_identical(colnames(fit$mat),
					 c("p.1", "q.2", "p.p.3", "q.q.4", "p.q.5"))
	expect_identical(rownames(fit$bal.tab),
					 c("p", "q", "p.p", "q.q", "p.q"))
	expect_identical(names(fit$grouping), c("linear", "squared", "two-way"))
	expect_equal(unname(fit$mat[, 3]), p^2, tolerance = 1e-8)
	expect_false(any(grepl("^X_", c(colnames(fit$mat), rownames(fit$bal.tab)))))
})

test_that("the treatment and outcome may be named like the generated covariates", {
	set.seed(51)
	n <- 400
	u <- rnorm(n); v <- rnorm(n)
	X1 <- rbinom(n, 1, plogis(0.6 * u))
	X2 <- 0.5 * X1 + u + v + rnorm(n)
	dat <- data.frame(u = u, v = v, X1 = X1, X2 = X2)

	fit <- hbal(Treat = "X1", X = c("u", "v"), Y = "X2", data = dat,
				print.level = -1)
	w <- fit$weights.co
	expect_lt(abs(mean(u[X1 == 1]) - sum(w * u[X1 == 0]) / sum(w)), 1e-3)
	expect_lt(abs(mean(v[X1 == 1]) - sum(w * v[X1 == 0]) / sum(w)), 1e-3)
	expect_identical(rownames(fit$bal.tab), c("u", "v"))
})

test_that("a covariate name that is not a column still gives the documented error", {
	set.seed(71)
	n <- 200
	a <- rnorm(n)
	treat <- rbinom(n, 1, plogis(0.5 * a))
	Y <- 0.5 * treat + a + rnorm(n)
	dat <- data.frame(a = a, treat = treat, Y = Y)

	expect_error(hbal(Treat = "treat", X = c("a", "NOPE"), Y = "Y", data = dat,
					  print.level = -1),
				 "not in the data")
})
