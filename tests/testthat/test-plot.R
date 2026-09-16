test_that("plot.hbal() draws the weight histogram with no deprecation warning", {
	old <- options(lifecycle_verbosity = "warning")
	on.exit(options(old), add = TRUE)

	set.seed(1984)
	N <- 200
	dat <- data.frame(treat = rbinom(N, 1, .5), X1 = rnorm(N), X2 = rbinom(N, 1, .5))
	dat$Y <- 0.5 * dat$treat + dat$X1 + dat$X2 + rnorm(N)
	out <- hbal(Treat = "treat", X = c("X1", "X2"), Y = "Y", data = dat, print.level = -1)

	# class= is required: a bare expect_no_warning() sets ignore_deprecation = TRUE
	# and would pass even while aes_string() is still warning.
	expect_no_warning(
		expect_message(p <- plot(out, type = "weight"), "sum\\(weights\\) normalized"),
		class = "lifecycle_warning_deprecated"
	)
	expect_s3_class(p, "ggplot")
	expect_identical(unname(vapply(p$layers, function(l) class(l$geom)[1], "")),
		c("GeomBar", "GeomDensity"))

	b <- ggplot2::ggplot_build(p)
	expect_equal(sum(b$data[[1]]$count), length(out$weights.co))
	rng <- range(log(out$weights.co))
	expect_lte(min(b$data[[1]]$xmin), rng[1])
	expect_gte(max(b$data[[1]]$xmax), rng[2])
})

test_that("plot.hbal() draws the covariate-balance panels with no deprecation warning", {
	old <- options(lifecycle_verbosity = "warning")
	on.exit(options(old), add = TRUE)
	pdf(NULL)
	on.exit(dev.off(), add = TRUE)

	set.seed(1984)
	N <- 200
	dat <- data.frame(treat = rbinom(N, 1, .5), X1 = rnorm(N), X2 = rbinom(N, 1, .5))
	dat$Y <- 0.5 * dat$treat + dat$X1 + dat$X2 + rnorm(N)
	out <- hbal(Treat = "treat", X = c("X1", "X2"), Y = "Y", data = dat, print.level = -1)

	expect_no_warning(g <- plot(out), class = "lifecycle_warning_deprecated")
	expect_s3_class(g, "gtable")
	expect_length(g$grobs, 2L)
	expect_length(g$grobs[[1]]$grobs, length(out$grouping))
})

test_that("plot.hbal() no longer references the deprecated aes_string()", {
	src <- deparse(getS3method("plot", "hbal"))
	expect_false(any(grepl("aes_string", src, fixed = TRUE)))
	expect_true(any(grepl('.data[["x"]]', src, fixed = TRUE)))
	expect_true(any(grepl('.data[["val"]]', src, fixed = TRUE)))
	expect_true(any(grepl('.data[["group"]]', src, fixed = TRUE)))
})
