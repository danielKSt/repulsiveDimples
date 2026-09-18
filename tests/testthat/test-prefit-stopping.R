# These run the whole estimation procedure, so they simulate. The fixture is the same
# deliberately tiny one the trust region loop tests use -- a 1x1 window, six patterns, a
# coarse r grid -- because what is being checked is when the pre-fitting stops and whether
# the bookkeeping around that is right, not statistical quality.

fixture <- function() {
  r <- seq(0, 0.3, by = 0.05)
  list(params = log(c(3, 0.15, 4)),
       repRange = 0.02,
       rho_hat = 12,
       K_hat = data.frame(r = r, border = pi * r^2),
       xlims = c(0, 1), ylims = c(0, 1), nSims = 6)
}

run_fit <- function(..., max.iter.prefit = 2, seed = 1) {
  f <- fixture()
  old <- options(mc.cores = 1)
  on.exit(options(old), add = TRUE)
  set.seed(seed)
  min_contrast_trust_region(params = f$params, parFreeIndex = c(1, 2, 3),
                            repRange = f$repRange, rho_hat = f$rho_hat, K_hat = f$K_hat,
                            xlims = f$xlims, ylims = f$ylims, nSims = f$nSims,
                            deltaInit = 0.2, eta = 0.05, deltaMax = 1, wq = c(1, 1/4),
                            subsection_count.main = 5, subsection_count.prefit = 5,
                            line_iterations.main = 1, line_iterations.prefit = 1,
                            eta_trust = 0.5, eta_converged = 0.99,
                            max.iter.prefit = max.iter.prefit, max.iter = 2, ...)
}

# Parameters at the end of a prefit cycle: the last block in it that actually ran.
cycle_end <- function(loop) {
  blocks <- Filter(Negate(is.null), loop)
  blocks[[length(blocks)]]$params
}

test_that("tolPrefitLoops = 0 runs every cycle maxPrefitLoops allows", {
  res <- run_fit(tolPrefitLoops = 0, maxPrefitLoops = 3)
  expect_length(res$prefitResults, 3L)
  expect_length(res$prefitMoves, 3L)
  expect_true(all(is.finite(res$prefitMoves)))
})

test_that("a tolerance nothing can beat stops after the first cycle", {
  # 10^6 is far above any movement a cycle on this fixture can produce, so the check has
  # to fire the first time it is reached.
  res <- run_fit(tolPrefitLoops = 10^6, maxPrefitLoops = 5)
  expect_length(res$prefitResults, 1L)
  expect_length(res$prefitMoves, 1L)
  # Trimmed, not padded: a caller reading length(prefitMoves) as the cycle count must not
  # be handed trailing NAs or NULLs.
  expect_false(anyNA(res$prefitMoves))
  expect_false(any(vapply(res$prefitResults, is.null, logical(1))))
})

test_that("prefitMoves is the distance the cycle actually moved the parameters", {
  f <- fixture()
  res <- run_fit(tolPrefitLoops = 0, maxPrefitLoops = 3)
  ends <- lapply(res$prefitResults, cycle_end)
  prev <- f$params
  for (k in seq_along(ends)) {
    expect_equal(res$prefitMoves[k], sqrt(sum((ends[[k]] - prev)^2)))
    prev <- ends[[k]]
  }
})

test_that("the loop stops on the first cycle that falls under the tolerance", {
  # Every cycle but the last must have moved at least the tolerance, or the loop ran on
  # past a cycle that should have ended it.
  res <- run_fit(tolPrefitLoops = 0.05, maxPrefitLoops = 6)
  moves <- res$prefitMoves
  if (length(moves) > 1) {
    expect_true(all(moves[-length(moves)] >= 0.05))
  }
  # The last entry says which rule stopped it: under the tolerance, or the cap binding.
  expect_true(moves[length(moves)] < 0.05 || length(moves) == 6L)
})

test_that("skipping the pre-fitting leaves both prefit slots NULL", {
  res <- run_fit(max.iter.prefit = 0)
  expect_null(res$prefitResults)
  expect_null(res$prefitMoves)

  res <- run_fit(maxPrefitLoops = 0)
  expect_null(res$prefitResults)
  expect_null(res$prefitMoves)
})

test_that("a negative tolerance is rejected rather than silently never firing", {
  expect_error(run_fit(tolPrefitLoops = -1), "tolPrefitLoops")
})

test_that("the method arguments reach both phases independently", {
  res <- run_fit(tolPrefitLoops = 0, maxPrefitLoops = 2,
                 method.prefit = "quadratic", method.main = "conjugate",
                 bisection_iterations.prefit = 4)
  expect_length(res$prefitResults, 2L)
  expect_true(all(is.finite(res$f_vals)))

  res <- run_fit(tolPrefitLoops = 0, maxPrefitLoops = 1,
                 method.prefit = "conjugate", method.main = "quadratic",
                 bisection_iterations.main = 4)
  expect_true(all(is.finite(res$f_vals)))
  expect_error(run_fit(method.main = "nonsense"))
})
