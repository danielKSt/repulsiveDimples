# trust_step and the line search below it are pure functions of the trust_function they
# are handed, so these tests supply a synthetic one and assert invariants: the step stays
# inside the trust region, never returns an iterate the importance sampling cannot
# support, and costs no more than the documented cap.

quad <- function(centre, ess_fun) {
  force(centre); force(ess_fun)
  function(x) list(f_est = sum((x - centre)^2), ess = ess_fun(x))
}

test_that("find_alpha_range treats delta as a radius, not a squared radius", {
  # Before 2026-09-10 this solved ||.||^2 = delta, so the region searched had radius
  # sqrt(delta) -- 2.2x too large at delta = 0.2, and the boundary tests elsewhere
  # compared plain distances against delta.
  for (delta in c(0.05, 0.2, 1, 2)) {
    r <- find_alpha_range(x_k = c(0, 0), x_0 = c(0, 0), p_k = c(1, 0), delta = delta)
    expect_equal(max(r), delta, tolerance = 1e-10)
    expect_equal(min(r), -delta, tolerance = 1e-10)
  }
})

test_that("find_alpha_range stays finite when the line is tangent to the boundary", {
  # The discriminant is 4a(delta^2 - d^2) and cannot truly be negative while x_k is
  # inside the region, but it can round below zero on the boundary. The NaN used to
  # reach seq() and abort the whole fit.
  r <- find_alpha_range(x_k = c(0.5, 0), x_0 = c(0, 0), p_k = c(0, 1),
                        delta = 0.5 * (1 - 1e-16))
  expect_true(all(is.finite(r)))
})

test_that("a zero-width range returns the starting point without evaluating anything", {
  calls <- 0
  tf <- function(x) { calls <<- calls + 1; list(f_est = sum(x^2), ess = 100) }
  res <- trust_line_steps(p_k = c(0, 1), x_k = c(0.5, 0), x_0 = c(0, 0), delta = 0.5,
                          trust_function = tf,
                          ess_bounds = repulsiveDimples:::ess_bounds_init(c(0, 0)),
                          eta_trust = 0, nSims = 100)
  expect_equal(res$x, c(0.5, 0))
  expect_equal(calls, 0)
})

test_that("the line search never returns a point failing the ESS requirement", {
  nSims <- 100
  # ESS collapses past a radius of 1, while the objective's minimum sits at 3.
  tf <- quad(c(3, 0), function(x) if (sqrt(sum(x^2)) <= 1) 0.9 * nSims else 0.2 * nSims)
  res <- trust_line_steps(p_k = c(1, 0), x_k = c(0, 0), x_0 = c(0, 0), delta = 5,
                          trust_function = tf,
                          ess_bounds = repulsiveDimples:::ess_bounds_init(c(0, 0)),
                          eta_trust = 0.6, nSims = nSims)
  expect_lte(sqrt(sum(res$x^2)), 1 + 1e-8)
  expect_gt(res$x[1], 0.9)          # still pushes up against the cliff
})

test_that("the line search falls back on x_k when the whole line is ruled out", {
  tf <- quad(c(0, 0), function(x) 10)     # every point fails at eta_trust = 0.6
  res <- trust_line_steps(p_k = c(1, 0), x_k = c(0.5, 0), x_0 = c(0, 0), delta = 5,
                          trust_function = tf,
                          ess_bounds = repulsiveDimples:::ess_bounds_init(c(0, 0)),
                          eta_trust = 0.6, nSims = 100)
  expect_equal(res$x, c(0.5, 0))
})

test_that("trust_step keeps its iterate inside the trust region", {
  set.seed(3)
  for (i in 1:25) {
    n <- sample(1:3, 1)
    centre <- rnorm(n, sd = 3)
    x0 <- rnorm(n, sd = 0.3)
    delta <- runif(1, 0.05, 1.5)
    r <- trust_step(x0, delta, quad(centre, function(x) 100), eta_trust = 0,
                    nSims = 100, max.directions = NULL)
    expect_lte(sqrt(sum((r$x_star - x0)^2)), delta * (1 + 1e-8))
  }
})

test_that("trust_step respects eta_trust at the iterate it returns", {
  nSims <- 100
  tf <- quad(c(3, 3), function(x) nSims * exp(-8 * abs(x[1]) - 0.2 * abs(x[2])))
  for (eta in c(0.3, 0.6, 0.8)) {
    r <- trust_step(c(0, 0), 0.5, tf, eta_trust = eta, nSims = nSims)
    expect_gt(r$ess_star / nSims, eta)
  }
})

test_that("subsection_count below 5 is rejected rather than indexing off the grid", {
  expect_error(trust_step(c(0, 0), 1, quad(c(1, 1), function(x) 100),
                          eta_trust = 0, nSims = 100, subsection_count = 4),
               "at least 5")
  expect_silent(trust_step(c(0, 0), 1, quad(c(1, 1), function(x) 100),
                           eta_trust = 0, nSims = 100, subsection_count = 5))
})

test_that("cost stays within the documented cap", {
  # (1 + max.directions*(n+1)) * (line_iterations+1) * subsection_count, plus the single
  # evaluation trust_step makes at the end to report f_pred and ess_star.
  for (cfg in list(c(11, 4, 5), c(11, 4, 1), c(7, 2, 2))) {
    n_eval <- 0
    tf <- function(x) { n_eval <<- n_eval + 1; list(f_est = sum(x^2), ess = 100) }
    trust_step(c(0, 0, 0), 1, tf, eta_trust = 0, nSims = 100,
               subsection_count = cfg[1], line_iterations = cfg[2], max.directions = cfg[3])
    expect_lte(n_eval, (1 + cfg[3] * 4) * (cfg[2] + 1) * cfg[1] + 1)
  }
})

test_that("a one-parameter step uses the finer grid the defaults have always paired with it", {
  n_eval <- 0
  tf <- function(x) { n_eval <<- n_eval + 1; list(f_est = (x - 0.3)^2, ess = 100) }
  trust_step(0, 1, tf, eta_trust = 0, nSims = 100,
             subsection_count = 11, line_iterations = 4)
  expect_equal(n_eval, 5 * 13 + 1)     # subsection_count + 2, over iterations + 1 rounds
})
