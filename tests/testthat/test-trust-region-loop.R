# These exercise the whole loop, so they simulate. The fixture is deliberately tiny -- a
# 1x1 window, six patterns, a coarse r grid -- because what is being checked is
# bookkeeping, not statistical quality.

fixture <- function() {
  r <- seq(0, 0.3, by = 0.05)
  list(params = log(c(3, 0.15, 4)),
       repRange = 0.02,
       rho_hat = 12,
       K_hat = data.frame(r = r, border = pi * r^2),
       xlims = c(0, 1), ylims = c(0, 1), nSims = 6)
}

run_loop <- function(..., seed = 1) {
  f <- fixture()
  old <- options(mc.cores = 1)
  on.exit(options(old), add = TRUE)
  set.seed(seed)
  trust_region_loop(params = f$params, parFreeIndex = c(1, 2, 3), repRange = f$repRange,
                    rho_hat = f$rho_hat, K_hat = f$K_hat, xlims = f$xlims, ylims = f$ylims,
                    nSims = f$nSims, deltaInit = 0.2, eta = 0.05, deltaMax = 1,
                    wq = c(1, 1/4), subsection_count = 5, line_iterations = 1, ...)
}

test_that("the iterate sequence and the contrast series stay in step", {
  res <- run_loop(eta_trust = 0.5, eta_converged = 0.99, max.iter = 4)
  expect_equal(nrow(res$x_seq), length(res$f_vals))
  expect_equal(ncol(res$x_seq), 3L)
  expect_true(all(is.finite(res$f_vals)))
  expect_equal(as.numeric(res$params), as.numeric(res$x_seq[nrow(res$x_seq), ]))
})

test_that("converging on ESS keeps the iterate it converged on", {
  # The loop used to break out before assigning update_eval$x, so a fit that converged
  # this way returned the parameters it started from and reported success. eta_converged
  # of 0 fires on the first step, which is exactly the path that was broken.
  f <- fixture()
  res <- run_loop(eta_trust = 0, eta_converged = 0, max.iter = 5)
  expect_equal(nrow(res$x_seq), 2L)
  expect_equal(length(res$f_vals), 2L)
  expect_equal(as.numeric(res$params), as.numeric(res$x_seq[2, ]))
  expect_true(all(is.finite(res$f_vals)))
})

test_that("validate.converged sends the converging iterate through the usual test", {
  # With validation on, the step is re-simulated and accepted or rejected like any other,
  # so the loop must not short-circuit and the bookkeeping must still line up.
  res <- run_loop(eta_trust = 0, eta_converged = 0, validate.converged = TRUE, max.iter = 5)
  expect_equal(nrow(res$x_seq), length(res$f_vals))
  expect_equal(as.numeric(res$params), as.numeric(res$x_seq[nrow(res$x_seq), ]))
  expect_length(res$update_evals, nrow(res$x_seq) - 1L)
})

test_that("fixed parameters are left alone", {
  f <- fixture()
  old <- options(mc.cores = 1)
  on.exit(options(old), add = TRUE)
  set.seed(2)
  res <- trust_region_loop(params = f$params, parFreeIndex = 2, repRange = f$repRange,
                           rho_hat = f$rho_hat, K_hat = f$K_hat, xlims = f$xlims,
                           ylims = f$ylims, nSims = f$nSims, deltaInit = 0.2, eta = 0.05,
                           deltaMax = 1, wq = c(1, 1/4), subsection_count = 5,
                           line_iterations = 1, eta_trust = 0.5, eta_converged = 0.99,
                           max.iter = 3)
  expect_equal(res$params[c(1, 3)], f$params[c(1, 3)])
  expect_equal(ncol(res$x_seq), 1L)
})

test_that("carry.directions runs end to end and does not disturb the bookkeeping", {
  res <- run_loop(eta_trust = 0.5, eta_converged = 0.99, max.iter = 4,
                  max.directions = NULL, carry.directions = TRUE)
  expect_equal(nrow(res$x_seq), length(res$f_vals))
  expect_true(all(is.finite(res$f_vals)))
})

test_that("method = quadratic runs the whole loop and keeps the bookkeeping", {
  res <- run_loop(eta_trust = 0.5, eta_converged = 0.99, max.iter = 4,
                  method = "quadratic", bisection_iterations = 4)
  expect_equal(nrow(res$x_seq), length(res$f_vals))
  expect_equal(ncol(res$x_seq), 3L)
  expect_true(all(is.finite(res$f_vals)))
  expect_equal(as.numeric(res$params), as.numeric(res$x_seq[nrow(res$x_seq), ]))
})

test_that("the quadratic method works on a block with one free parameter", {
  # The mu prefit block has a single free parameter, where a quadratic needs three points
  # and there are no cross terms to place.
  f <- fixture()
  old <- options(mc.cores = 1)
  on.exit(options(old), add = TRUE)
  set.seed(3)
  res <- trust_region_loop(params = f$params, parFreeIndex = 3, repRange = f$repRange,
                           rho_hat = f$rho_hat, K_hat = f$K_hat, xlims = f$xlims,
                           ylims = f$ylims, nSims = f$nSims, deltaInit = 0.2, eta = 0.05,
                           deltaMax = 1, wq = c(1, 1/4), eta_trust = 0.5,
                           eta_converged = 0.99, max.iter = 3, method = "quadratic",
                           bisection_iterations = 4)
  expect_equal(res$params[c(1, 2)], f$params[c(1, 2)])
  expect_equal(ncol(res$x_seq), 1L)
  expect_true(all(is.finite(res$f_vals)))
})

test_that("carry.ess_bounds runs end to end without disturbing the bookkeeping", {
  for (meth in c("conjugate", "quadratic")) {
    res <- run_loop(eta_trust = 0.5, eta_converged = 0.99, max.iter = 4, method = meth,
                    bisection_iterations = 4, carry.ess_bounds = TRUE)
    expect_equal(nrow(res$x_seq), length(res$f_vals))
    expect_true(all(is.finite(res$f_vals)))
  }
})
