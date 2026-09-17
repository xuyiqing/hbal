# Tests for balanceData(), the accessor offered in GitHub PR #4
# (joshuafayallen). It must return exactly the data plot.hbal() draws, so a
# caller can rebuild the balance figure with ggplot2.

test_that("balanceData() returns the documented shape", {
	set.seed(1984)
	N <- 300
	X1 <- rnorm(N)
	X2 <- rbinom(N, size = 1, prob = .5)
	treat <- rbinom(N, 1, prob = 0.5)
	Y <- 0.5 * treat + X1 + X2 + rnorm(N)
	dat <- data.frame(treat = treat, X1 = X1, X2 = X2, Y = Y)
	out <- hbal(Treat = "treat", X = c("X1", "X2"), Y = "Y", data = dat,
				print.level = -1)

	b <- balanceData(out)
	expect_s3_class(b, "data.frame")
	expect_named(b, c("term", "covar.group", "adjustment", "std.diff",
					  "tr.mean", "co.mean", "w.co.mean"))
	expect_equal(nrow(b), 2L * nrow(out$bal.tab))
	expect_type(b$term, "character")
	expect_s3_class(b$covar.group, "factor")
	expect_identical(levels(b$adjustment), c("before", "after"))
	expect_identical(rownames(b), as.character(seq_len(nrow(b))))
})

test_that("balanceData() matches what plot.hbal() draws", {
	set.seed(202)
	n <- 400
	aa <- rnorm(n)
	bb <- rnorm(n)
	treat <- rbinom(n, 1, plogis(0.6 * aa))
	Y <- 0.5 * treat + aa + bb + rnorm(n)
	dat <- data.frame(aa = aa, bb = bb, treat = treat, Y = Y)
	out <- hbal(Treat = "treat", X = c("aa", "bb"), Y = "Y", data = dat,
				expand.degree = 2, print.level = -1)

	b <- balanceData(out)
	groups <- names(out$grouping)
	expect_identical(b$term, rep(colnames(out$mat), 2))
	expect_identical(as.character(b$covar.group),
					 rep(rep(groups, out$grouping), 2))
	expect_equal(unname(b$std.diff),
				 unname(c(out$bal.tab[, 4], out$bal.tab[, 5])))
	expect_identical(levels(b$covar.group), groups)
	expect_equal(unname(b$tr.mean), unname(rep(out$bal.tab[, 1], 2)))
})

test_that("balanceData() rejects anything that is not an hbal object", {
	expect_error(balanceData(list()), "must be an hbal object")
	expect_error(balanceData(data.frame(a = 1)), "must be an hbal object")
})
