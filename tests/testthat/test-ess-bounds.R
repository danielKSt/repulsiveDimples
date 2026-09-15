# The ess_bounds box screens out evaluations the importance sampling cannot support.
# The tests below are about its invariants rather than particular numbers, because the
# one bug this machinery had was an invariant violation: recording a failure against
# every parameter the failing point moved ratcheted the box shut around the search, to
# the point where a returned iterate sat outside its own bounds.

init     <- repulsiveDimples:::ess_bounds_init
record   <- repulsiveDimples:::ess_bounds_record
feasible <- repulsiveDimples:::ess_bounds_feasible

test_that("a fresh box rules nothing out", {
  b <- init(c(0, 0, 0))
  expect_true(all(is.infinite(b$lower) & b$lower < 0))
  expect_true(all(is.infinite(b$upper) & b$upper > 0))
  expect_equal(b$pass_lower, c(0, 0, 0))
  expect_equal(b$pass_upper, c(0, 0, 0))
  expect_true(feasible(b, c(1e6, -1e6, 3)))
})

test_that("a failure bounds only the parameters that moved, on the side they moved to", {
  b <- record(init(c(0, 0, 0)), c(2, 0, -3), passed = FALSE)
  expect_equal(b$upper, c(2, Inf, Inf))
  expect_equal(b$lower, c(-Inf, -Inf, -3))
})

test_that("bounds only ever tighten", {
  b <- record(init(c(0, 0, 0)), c(2, 0, 0), passed = FALSE)
  b <- record(b, c(5, 0, 0), passed = FALSE)      # looser: must not widen
  expect_equal(b$upper[1], 2)
  b <- record(b, c(1, 0, 0), passed = FALSE)      # tighter: must take
  expect_equal(b$upper[1], 1)
})

test_that("a bound is exclusive: a point sitting on it is ruled out", {
  b <- record(init(c(0, 0, 0)), c(2, 0, 0), passed = FALSE)
  expect_false(feasible(b, c(2, 0, 0)))
  expect_false(feasible(b, c(3, 0, 0)))
  expect_true(feasible(b, c(1.999, 0, 0)))
})

test_that("a value that has passed is never bounded away (the ratchet guard)", {
  # A breakdown is a property of the whole vector. Blaming every parameter the failing
  # point moved lets a failure driven by parameter 2 bound parameter 1 at a value that
  # already passed, which closes the box around the search and eventually excludes the
  # search's own iterate.
  b <- init(c(0, 0))
  b <- record(b, c(0.45, 0.0), passed = TRUE)
  b <- record(b, c(0.45, 0.9), passed = FALSE)    # failure driven by parameter 2
  expect_equal(b$upper, c(Inf, 0.9))
  expect_true(feasible(b, c(0.45, 0.0)))
  b <- record(b, c(0.6, 0.0), passed = FALSE)     # a genuine parameter 1 failure
  expect_equal(b$upper, c(0.6, 0.9))
})

test_that("under the real evaluation protocol, a passed point stays feasible", {
  # trust_line_evaluate screens a point against the box BEFORE evaluating it, so a point
  # outside the box is never evaluated and so never recorded as passed. The invariant is
  # therefore about points admitted under that protocol, and this test follows it: check
  # first, and only record what the screen would have let through. Skipping the check and
  # recording arbitrary points does violate the invariant, but that state is unreachable.
  set.seed(4)
  b <- init(c(0, 0, 0))
  passed <- list()
  admitted <- 0
  for (i in 1:400) {
    x <- rnorm(3, sd = 0.5)
    if (!feasible(b, x)) next
    admitted <- admitted + 1
    ok <- sum(abs(x)) < 0.8
    b <- record(b, x, passed = ok)
    if (ok) passed[[length(passed) + 1]] <- x
    for (p in passed) expect_true(feasible(b, p))
  }
  expect_gt(admitted, 20)
  expect_gt(length(passed), 0)
})

test_that("x_0 is always feasible, so a line search always has a fallback", {
  # trust_line_steps returns x_k when the whole line is ruled out, which is only sound
  # if the centre itself can never be bounded away.
  set.seed(9)
  b <- init(c(0, 0, 0))
  for (i in 1:200) {
    x <- rnorm(3, sd = 0.5)
    if (!feasible(b, x)) next
    b <- record(b, x, passed = sum(abs(x)) < 0.8)
    expect_true(feasible(b, c(0, 0, 0)))
  }
})
