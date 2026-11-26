# estimate_K_lambda_baseline finds the pairs within the largest radius with closepairs.
# The reference here goes through every ordered pair of points directly, so it shares
# nothing with the implementation. The window is not square, so that a mix-up of the x
# and y sides of the translation correction shows.

K_reference <- function(x, y, xlim, ylim, r_vec) {
  s <- numeric(length(r_vec))
  for (i in seq_along(x)) for (j in seq_along(x)) if (i != j) {
    dx <- x[i] - x[j]
    dy <- y[i] - y[j]
    w <- 1/((diff(xlim) - abs(dx))*(diff(ylim) - abs(dy)))
    s <- s + w*(sqrt(dx^2 + dy^2) <= r_vec)
  }
  s
}

as_pattern <- function(x, y, xlim, ylim) {
  list(thinned = data.frame(x = x, y = y), xlim_thinned = xlim, ylim_thinned = ylim)
}

test_that("the K estimate matches a direct sum over every pair", {
  xlim <- c(0, 4)
  ylim <- c(-1, 2)
  r_vec <- seq(from = 0, to = 1.5, by = 0.05)
  set.seed(5)
  for (n in c(2, 3, 40, 150)) {
    x <- stats::runif(n, xlim[1], xlim[2])
    y <- stats::runif(n, ylim[1], ylim[2])
    got <- estimate_K_lambda_baseline(as_pattern(x, y, xlim, ylim), r_vec)
    expect_equal(got$r, r_vec)
    expect_equal(got$trans, K_reference(x, y, xlim, ylim, r_vec), tolerance = 1e-12)
  }
})

test_that("a pair at exactly the largest radius is counted there", {
  xlim <- c(0, 4)
  ylim <- c(0, 3)
  r_vec <- c(0, 0.5, 1, 1.5)
  x <- c(0.5, 2.0, 3.5)
  y <- c(0.5, 0.5, 2.5)
  got <- estimate_K_lambda_baseline(as_pattern(x, y, xlim, ylim), r_vec)
  expect_gt(got$trans[4], 0)
  expect_equal(got$trans, K_reference(x, y, xlim, ylim, r_vec), tolerance = 1e-12)
})

test_that("patterns without close pairs, or with fewer than two points, give zero", {
  r_vec <- seq(from = 0, to = 1, by = 0.25)
  far <- as_pattern(c(0.5, 3.5), c(0.5, 2.5), c(0, 4), c(0, 3))
  expect_equal(estimate_K_lambda_baseline(far, r_vec)$trans, rep(0, length(r_vec)))
  one <- as_pattern(1, 1, c(0, 4), c(0, 3))
  expect_equal(estimate_K_lambda_baseline(one, r_vec)$trans, rep(0, length(r_vec)))
  none <- list(thinned = data.frame(), xlim_thinned = c(0, 4), ylim_thinned = c(0, 3))
  expect_equal(estimate_K_lambda_baseline(none, r_vec)$trans, rep(0, length(r_vec)))
})
