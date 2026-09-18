test_that("att() returns the documented output shape with no warnings or messages", {
	set.seed(1984)
	N <- 100
	X1 <- rnorm(N)
	X2 <- rbinom(N, size = 1, prob = .5)
	treat <- rbinom(N, 1, prob = 0.5)
	y <- 0.5 * treat + X1 + X2 + rnorm(N)
	dat <- data.frame(treat = treat, X1 = X1, X2 = X2, Y = y)

	out <- hbal(Treat = "treat", X = c("X1", "X2"), Y = "Y", data = dat)

	expect_no_warning(res <- att(out))
	expect_no_message(att(out))
	expect_s3_class(res, "data.frame")
	expect_false(inherits(res, "tbl_df"))
	expect_equal(dim(res), c(1L, 7L))
	expect_named(res, c("Estimate", "Std. Error", "t value", "Pr(>|t|)",
						 "CI Lower", "CI Upper", "DF"))
	expect_identical(rownames(res), "treat")

	expect_no_warning(att(out, method = "lm_lin"))
	expect_no_warning(att(out, dr = FALSE))
	expect_no_warning(att(out, method = "lm_robust", se_type = "HC2"))
})

test_that(".att_tidy_select() handles a tibble-classed input without warning", {
	stub <- data.frame(
		term = c("(Intercept)", "treat", "X1"),
		estimate = c(1.1, 0.5, 0.2),
		std.error = c(0.1, 0.05, 0.02),
		statistic = c(11, 10, 10),
		p.value = c(0.001, 0.001, 0.001),
		conf.low = c(0.9, 0.4, 0.16),
		conf.high = c(1.3, 0.6, 0.24),
		df = c(97, 97, 97),
		outcome = c("Y", "Y", "Y"),
		stringsAsFactors = FALSE
	)
	class(stub) <- c("tbl_df", "tbl", "data.frame")

	expect_no_warning(res <- hbal:::.att_tidy_select(stub, "treat"))
	expect_s3_class(res, "data.frame")
	expect_false(inherits(res, "tbl_df"))
	expect_equal(dim(res), c(1L, 7L))
	expect_named(res, c("Estimate", "Std. Error", "t value", "Pr(>|t|)",
						 "CI Lower", "CI Upper", "DF"))
	expect_identical(rownames(res), "treat")
	expect_equal(as.numeric(res),
				 c(0.5, 0.05, 10, 0.001, 0.4, 0.6, 97))
})

test_that(".att_tidy_select() errors loudly on a missing or duplicated term", {
	stub <- data.frame(
		term = c("(Intercept)", "treat", "X1"),
		estimate = c(1.1, 0.5, 0.2), std.error = c(0.1, 0.05, 0.02),
		statistic = c(11, 10, 10), p.value = c(0.001, 0.001, 0.001),
		conf.low = c(0.9, 0.4, 0.16), conf.high = c(1.3, 0.6, 0.24),
		df = c(97, 97, 97), outcome = c("Y", "Y", "Y"),
		stringsAsFactors = FALSE
	)

	expect_error(hbal:::.att_tidy_select(stub, "nonexistent"))

	dup <- rbind(stub, stub[2, ])
	expect_error(hbal:::.att_tidy_select(dup, "treat"))

	missing_col <- stub[, setdiff(colnames(stub), "conf.low")]
	expect_error(hbal:::.att_tidy_select(missing_col, "treat"))
})
