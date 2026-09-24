# The conjugate direction set is carried across trust region iterations under
# carry.directions. Restarting it from the coordinate axes every step discards Powell's
# accumulation and makes every step open along free parameter 1.

dsi <- repulsiveDimples:::direction_set_init

# Direction sets only exist in the conjugate step, so every trust_step call below asks for
# it by name. The default has been the quadratic step since 2026-09-24, and without the
# pin these would silently run a method that returns no direction set at all.

test_that("direction_set_init falls back to the axes on anything unusable", {
  expect_equal(dsi(NULL, 3), diag(3))
  expect_equal(dsi(diag(2), 3), diag(3))                                   # wrong shape
  expect_equal(dsi(matrix(c(1, 0, 2, 0, 0, 0, 0, 1, 0), 3, 3), 3), diag(3))  # singular
  expect_equal(dsi(matrix(c(NA, rep(0, 8)), 3, 3), 3), diag(3))            # non-finite
  expect_equal(dsi(matrix(0, 3, 3), 3), diag(3))
})

test_that("direction_set_init passes a usable set straight through", {
  set.seed(2)
  d <- qr.Q(qr(matrix(rnorm(9), 3, 3)))
  expect_equal(dsi(d, 3), d)
})

test_that("trust_step returns a direction set its own guard accepts", {
  set.seed(5)
  tf <- function(x) list(f_est = as.numeric(crossprod(x - c(2, 1)) +
                                             0.5 * (x[1] - 2) * (x[2] - 1)), ess = 100)
  r <- trust_step(c(0, 0), 0.5, tf, eta_trust = 0, nSims = 100, method = "conjugate", max.directions = NULL)
  expect_true(is.matrix(r$directions))
  expect_equal(dim(r$directions), c(2L, 2L))
  expect_equal(dsi(r$directions, 2), r$directions)
})

test_that("a carried set is actually used, and an uncarried one is not", {
  # Feeding a rotated set in must change the first direction searched, and so the path.
  tf <- function(x) list(f_est = as.numeric(crossprod(x - c(2, 1))), ess = 100)
  rot <- matrix(c(0, 1, 1, 0), 2, 2)             # swaps which parameter is searched first
  a <- trust_step(c(0, 0), 0.3, tf, eta_trust = 0, nSims = 100, method = "conjugate",
                  max.directions = 1, directions = NULL)
  b <- trust_step(c(0, 0), 0.3, tf, eta_trust = 0, nSims = 100, method = "conjugate",
                  max.directions = 1, directions = rot)
  expect_false(isTRUE(all.equal(a$x_star, b$x_star)))
})

test_that("the stored direction is a unit vector", {
  # Not cosmetic: unnormalised, the stored net displacements shrink as a fit converges
  # until det() underflows to zero and the set silently resets to the axes, which would
  # make carrying quietly equivalent to not carrying.
  set.seed(8)
  tf <- function(x) list(f_est = as.numeric(crossprod(x - c(3, 2)) +
                                             0.7 * (x[1] - 3) * (x[2] - 2)), ess = 100)
  dirs <- NULL
  x <- c(0, 0)
  for (i in 1:5) {
    r <- trust_step(x, 0.3, tf, eta_trust = 0, nSims = 100, method = "conjugate",
                    max.directions = NULL, directions = dirs)
    x <- r$x_star
    dirs <- r$directions
    expect_gt(abs(det(dirs)), 1e-8)
    norms <- sqrt(colSums(dirs^2))
    expect_true(all(abs(norms - 1) < 1e-8 | abs(norms) < 1e-12))
  }
})

test_that("carrying reaches a correlated quadratic's minimum no worse, and no dearer", {
  A <- matrix(c(4, 3, 0, 3, 4, 0, 0, 0, 2), 3, 3)
  bvec <- c(-1, 2, 0.5)
  fq <- function(x) as.numeric(0.5 * t(x) %*% A %*% x + sum(bvec * x))
  opt <- solve(A, -bvec)
  run <- function(carry) {
    n <- 0
    x <- c(0, 0, 0); dirs <- NULL
    for (i in 1:6) {
      r <- trust_step(x, 0.4, function(z) { n <<- n + 1; list(f_est = fq(z), ess = 100) },
                      eta_trust = 0, nSims = 100, method = "conjugate", max.directions = NULL,
                      directions = if (carry) dirs else NULL)
      x <- r$x_star
      if (carry) dirs <- r$directions
    }
    list(err = sqrt(sum((x - opt)^2)), n = n)
  }
  off <- run(FALSE); on <- run(TRUE)
  expect_lte(on$err, off$err + 1e-8)
  expect_lte(on$n, off$n)
})
