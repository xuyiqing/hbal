# Shared fixtures for tests/testthat/test-att-aipw.R (auto-sourced by testthat).
# print.level = -1 only silences hbal()'s informational message()s; it changes
# no computed value.

make_toy <- function(seed, n1 = 150, n0 = 150, expand.degree = 1, cv = FALSE) {
  set.seed(seed)  # test code — set.seed() here is fine, this is not package code
  n <- n1 + n0
  X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rbinom(n, 1, 0.5)
  Tr <- c(rep(1, n1), rep(0, n0))
  Y <- 1 + 0.5 * Tr + X1 + 0.8 * X2 - 0.6 * X3 + rnorm(n)
  dat <- data.frame(Tr = Tr, X1 = X1, X2 = X2, X3 = X3, Y = Y)
  hbal::hbal(Treat = "Tr", X = c("X1", "X2", "X3"), Y = "Y", data = dat,
             expand.degree = expand.degree, cv = cv, print.level = -1)
}

# Small fixture: n0 = 10 < 4p = 12 forces the K = 1 fallback (no cross-fitting),
# so the closed-form agreement under exact balance applies without cross-fitting
# slack. Used by the closed-form test and the fold-fallback test.
make_small_toy <- function() {
  set.seed(777)
  n1 <- 20; n0 <- 10
  X1 <- rnorm(n1 + n0); X2 <- rnorm(n1 + n0); X3 <- rbinom(n1 + n0, 1, 0.5)
  Tr <- c(rep(1, n1), rep(0, n0))
  Y  <- 1 + 0.5 * Tr + X1 + 0.8 * X2 - 0.6 * X3 + rnorm(n1 + n0)
  dat <- data.frame(Tr = Tr, X1 = X1, X2 = X2, X3 = X3, Y = Y)
  hbal::hbal(Treat = "Tr", X = c("X1", "X2", "X3"), Y = "Y", data = dat,
             expand.degree = 1, print.level = -1)
}
