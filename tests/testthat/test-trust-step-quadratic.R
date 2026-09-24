# The model fitting and minimisation are checked against quadratics whose answer is known
# in closed form, so a failure points at the algebra rather than at the optimizer. The ESS
# machinery is checked against a synthetic trust function with a feasible region set by
# hand, which keeps these fast and makes the expected box exact.

test_that("the interpolation recovers a known quadratic exactly", {
  set.seed(5)
  for (n in 1:3) {
    g <- round(rnorm(n), 2)
    H <- matrix(round(rnorm(n * n), 2), n, n)
    H <- H + t(H)
    fq <- function(s) 1.7 + sum(g * s) + 0.5 * sum(s * as.vector(H %*% s))

    box <- list(lower = rep(-0.3, n), upper = rep(0.25, n))
    S <- matrix(0, nrow = 1, ncol = n)
    for (i in 1:n) for (side in c(box$upper[i], box$lower[i])) {
      s <- rep(0, n); s[i] <- side; S <- rbind(S, s)
    }
    if (n > 1) for (i in 1:(n - 1)) for (j in (i + 1):n) {
      s <- rep(0, n); s[i] <- box$upper[i]; s[j] <- box$upper[j]; S <- rbind(S, s)
    }
    rownames(S) <- NULL

    expect_equal(nrow(S), (n + 1) * (n + 2) / 2)
    m <- quadratic_model_fit(S = S, f = apply(S, 1, fq), n = n)
    expect_false(is.null(m))
    expect_equal(m$g, as.numeric(g))
    expect_equal(m$H, H)
  }
})

test_that("a set that cannot determine a quadratic is refused rather than fitted", {
  # Too few usable points, and a set that has collapsed onto a plane, both leave the
  # interpolation system without a unique solution.
  S <- matrix(0, nrow = 10, ncol = 3)
  expect_null(quadratic_model_fit(S = S[1:4, , drop = FALSE], f = rep(1, 4), n = 3))
  expect_null(quadratic_model_fit(S = S, f = rep(1, 10), n = 3))

  S2 <- rbind(c(0,0,0), c(1,0,0), c(-1,0,0), c(0,1,0), c(0,-1,0),
              c(0,0,0), c(0,0,0), c(1,1,0), c(0,0,0), c(0,0,0))
  expect_null(quadratic_model_fit(S = S2, f = seq_len(10), n = 3))
})

test_that("the model minimiser respects the box and the ball", {
  # A linear model pushes straight to the boundary, so the active constraint is whichever
  # of the two is tighter, and the answer is known.
  m <- list(g = c(-1, 0), H = matrix(0, 2, 2))
  s <- quadratic_model_minimise(m, lower = c(-5, -5), upper = c(0.3, 5), delta = 10,
                                starts = rbind(c(0, 0), c(1, 1)))
  expect_equal(s[1], 0.3)

  s <- quadratic_model_minimise(m, lower = c(-5, -5), upper = c(5, 5), delta = 0.5,
                                starts = rbind(c(0, 0), c(1, 1)))
  expect_lte(sqrt(sum(s^2)), 0.5 + 1e-9)
  expect_equal(sqrt(sum(s^2)), 0.5)

  # An indefinite model must still come back inside the region rather than running away.
  m2 <- list(g = c(0, 0), H = diag(c(-2, -2)))
  s <- quadratic_model_minimise(m2, lower = c(-1, -1), upper = c(1, 1), delta = 0.4,
                                starts = rbind(c(0, 0), c(0.1, 0.1)))
  expect_lte(sqrt(sum(s^2)), 0.4 + 1e-9)
})

# A trust function whose ESS is 1 inside a hand-set box and 0 outside it, so the bisection
# has an exact answer to find.
box_trust <- function(halfwidth, nSims = 100) {
  function(x) {
    inside <- all(abs(x) <= halfwidth + 1e-12)
    list(f_est = sum(x^2), ess = if (inside) nSims else 0)
  }
}

test_that("the ray search finds the crossing, and stops at delta when delta binds", {
  tf <- box_trust(c(0.05, 0.4, 0.4))
  st <- list(ess_bounds = ess_bounds_init(rep(0, 3)), evals = list())

  # Coordinate 1 crosses inside the radius, so it is bisected to the halfwidth.
  r <- ess_ray_limit(rep(0, 3), c(1, 0, 0), t_hi = 0.2, trust_function = tf,
                     eta_trust = 0.5, nSims = 100, bisection_iterations = 20, state = st)
  expect_lte(r$t, 0.05)
  expect_gt(r$t, 0.05 - 1e-4)

  # Coordinate 2 does not cross inside the radius, so the radius is returned and the
  # search costs a single evaluation.
  n_called <- 0
  counting <- function(x) { n_called <<- n_called + 1; tf(x) }
  r <- ess_ray_limit(rep(0, 3), c(0, 1, 0), t_hi = 0.2, trust_function = counting,
                     eta_trust = 0.5, nSims = 100, bisection_iterations = 20, state = st)
  expect_equal(r$t, 0.2)
  expect_equal(n_called, 1L)
})

test_that("the box is the per-coordinate crossing, capped at the trust radius", {
  tf <- box_trust(c(0.05, 0.4, 0.12))
  st <- list(ess_bounds = ess_bounds_init(rep(0, 3)), evals = list())
  b <- ess_box_bisect(rep(0, 3), delta = 0.2, trust_function = tf, eta_trust = 0.5,
                      nSims = 100, bisection_iterations = 20, state = st)
  expect_equal(b$box$upper, c(0.05, 0.2, 0.12), tolerance = 1e-4)
  expect_equal(b$box$lower, c(-0.05, -0.2, -0.12), tolerance = 1e-4)
  expect_true(all(b$box$lower <= 0) && all(b$box$upper >= 0))
})

test_that("every interpolation point ends up inside the supported region", {
  # The corner of the box is deliberately outside the feasible set here, so the repair
  # has to fire and must not leave a hole behind.
  halfwidth <- c(0.1, 0.1, 0.1)
  tf <- function(x) list(f_est = sum(x^2),
                         ess = if (sum(abs(x)) <= 0.12) 100 else 0)
  st <- list(ess_bounds = ess_bounds_init(rep(0, 3)), evals = list())
  b <- ess_box_bisect(rep(0, 3), delta = 0.5, trust_function = tf, eta_trust = 0.5,
                      nSims = 100, bisection_iterations = 25, state = st)
  set <- quadratic_interp_set(rep(0, 3), b$box, delta = 0.5, trust_function = tf,
                              eta_trust = 0.5, nSims = 100, bisection_iterations = 25,
                              state = b$state)
  expect_equal(nrow(set$S), 10L)
  expect_true(all(is.finite(set$f)))
  expect_true(all(apply(set$S, 1, function(s) tf(s)$ess) > 50))
})

test_that("the quadratic step returns the documented shape and a supported iterate", {
  tf <- function(x) {
    inside <- sqrt(sum(x^2)) <= 0.3
    list(f_est = sum((x - c(0.05, -0.02, 0.01))^2), ess = if (inside) 100 else 0)
  }
  r <- trust_step_quadratic(x_0 = rep(0, 3), delta = 0.2, trust_function = tf,
                            eta_trust = 0.5, nSims = 100)
  expect_named(r, c("x_star", "f_pred", "ess_star", "ess_bounds", "directions", "box"),
               ignore.order = TRUE)
  expect_null(r$directions)
  expect_length(r$x_star, 3L)
  expect_gt(r$ess_star / 100, 0.5)
  # The minimum of this trust function is inside the region, so the step must find it.
  expect_equal(r$x_star, c(0.05, -0.02, 0.01), tolerance = 1e-2)
  expect_lt(r$f_pred, tf(rep(0, 3))$f_est)
})

test_that("the step never returns a point worse than the ones it evaluated", {
  # A deliberately non-quadratic trust function, so the model is wrong and the fallback
  # on the best measured point is what has to hold.
  tf <- function(x) list(f_est = sum(abs(x)^(1/2)) + 10 * sin(6 * x[1]),
                         ess = if (sqrt(sum(x^2)) <= 0.25) 100 else 0)
  r <- trust_step_quadratic(x_0 = rep(0, 3), delta = 0.2, trust_function = tf,
                            eta_trust = 0.5, nSims = 100)
  expect_true(is.finite(r$f_pred))
  expect_lte(r$f_pred, tf(rep(0, 3))$f_est + 1e-9)
})

test_that("interp_fraction is validated", {
  tf <- function(x) list(f_est = sum(x^2), ess = 100)
  expect_error(trust_step_quadratic(rep(0, 3), 0.2, tf, 0.5, 100, interp_fraction = 0),
               "interp_fraction")
  expect_error(trust_step_quadratic(rep(0, 3), 0.2, tf, 0.5, 100, interp_fraction = 1.5),
               "interp_fraction")
})

test_that("trust_step dispatches on method and both methods return the same shape", {
  tf <- function(x) list(f_est = sum((x - c(0.04, 0.01, -0.02))^2),
                         ess = if (sqrt(sum(x^2)) <= 0.3) 100 else 0)
  a <- trust_step(rep(0, 3), 0.2, tf, 0.5, 100, method = "conjugate",
                  subsection_count = 5, line_iterations = 1)
  b <- trust_step(rep(0, 3), 0.2, tf, 0.5, 100, method = "quadratic")
  for (nm in c("x_star", "f_pred", "ess_star", "ess_bounds", "directions")) {
    expect_true(nm %in% names(a))
    expect_true(nm %in% names(b))
  }
  expect_error(trust_step(rep(0, 3), 0.2, tf, 0.5, 100, method = "nonsense"))
  # The default has been the quadratic step since 2026-09-24.
  expect_equal(trust_step(rep(0, 3), 0.2, tf, 0.5, 100)$x_star, b$x_star)
})

test_that("carried bounds are used at the same centre and discarded at any other", {
  x0 <- c(0.1, -0.2, 0.3)
  b <- ess_bounds_init(x0)
  b <- ess_bounds_record(b, x0 + c(0.5, 0, 0), passed = FALSE)
  expect_equal(ess_bounds_carried(b, x0)$upper[1], x0[1] + 0.5)

  # A different centre, the wrong length, and nothing at all all give a fresh set.
  expect_equal(ess_bounds_carried(b, x0 + 1)$upper, rep(Inf, 3))
  expect_equal(ess_bounds_carried(b, c(0, 0))$upper, rep(Inf, 2))
  expect_equal(ess_bounds_carried(NULL, x0)$upper, rep(Inf, 3))
  expect_equal(ess_bounds_carried(NULL, x0)$x_0, x0)
})

test_that("carrying bounds across a rejected step costs evaluations, not accuracy", {
  # The situation trust_region_loop creates when a step is rejected: same centre, same
  # trust function, a smaller radius. The iterate must not move, and the second step must
  # not pay again for the region the first one already mapped.
  halfwidth <- c(0.06, 0.4, 0.1)
  mk <- function() {
    n <- 0
    list(f = function(x) {
           n <<- n + 1
           list(f_est = sum((x - c(0.03, 0.02, -0.04))^2),
                ess = if (all(abs(x) <= halfwidth)) 100 else 0)
         },
         n = function() n)
  }
  for (meth in c("quadratic", "conjugate")) {
    w0 <- mk()
    first <- trust_step(rep(0, 3), 0.2, w0$f, 0.5, 100, method = meth)

    wA <- mk()
    fresh <- trust_step(rep(0, 3), 0.18, wA$f, 0.5, 100, method = meth, ess_bounds = NULL)
    wB <- mk()
    carried <- trust_step(rep(0, 3), 0.18, wB$f, 0.5, 100, method = meth,
                          ess_bounds = first$ess_bounds)

    expect_equal(fresh$x_star, carried$x_star)
    expect_lte(wB$n(), wA$n())
  }
})
