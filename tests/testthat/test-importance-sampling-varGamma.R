# At nu = 1/2 the Variance gamma kernel is the exponential kernel
# exp(-r/scale)/(2*pi*scale^2), which gives the kernel sums something exact to be checked
# against.

scale <- 0.3
xl <- c(0, 3)
yl <- c(0, 3)
k_exp <- function(r) exp(-r/scale)/(2*pi*scale^2)

no_daughters <- data.frame(x = numeric(0), y = numeric(0))

# Direct 2D integral of the exponential kernel over the window, with each axis split at
# the parent so that the kernel's peak sits on a breakpoint.
p_W_2d <- function(cx, cy) {
  split <- function(lim, at) sort(unique(c(lim, min(max(at, lim[1]), lim[2]))))
  piecewise <- function(f, br, tol) {
    sum(vapply(seq_len(length(br) - 1), function(i) {
      stats::integrate(f, br[i], br[i + 1], rel.tol = tol, abs.tol = 1e-14)$value
    }, numeric(1)))
  }
  inner <- function(y) vapply(y, function(yy) {
    piecewise(function(x) k_exp(sqrt((x - cx)^2 + (yy - cy)^2)), split(xl, cx), 1e-9)
  }, numeric(1))
  piecewise(inner, split(yl, cy), 1e-8)
}

patterns <- function(n = 4) {
  set.seed(3)
  lapply(seq_len(n), function(i) {
    rVarGamma_matern_thinned(kappa = 1, scale = scale, mu = 3, nu = 1/2, repulsionRange = 0.1,
                             xlims = xl, ylims = yl, saveparents = TRUE)
  })
}


test_that("the window probabilities match a direct integral of the kernel", {
  # a parent in the centre, just inside a corner, and just outside an edge
  parent <- data.frame(x = c(1.5, 0.05, -0.1), y = c(1.5, 0.05, 1.5))
  for (j in seq_len(nrow(parent))) {
    sums <- varGamma_daughter_kernel_sums(scale = scale, nu = 1/2, parent = parent[j, ],
                                          daughter = no_daughters, xlim = xl, ylim = yl)
    expect_equal(sums$S_pW, p_W_2d(parent$x[j], parent$y[j]), tolerance = 1e-9)
  }
})

test_that("the log kernel sums match the exponential kernel", {
  set.seed(1)
  parent <- data.frame(x = stats::runif(20, -1, 4), y = stats::runif(20, -1, 4))
  daughter <- data.frame(x = stats::runif(50, 0, 3), y = stats::runif(50, 0, 3))
  sums <- varGamma_daughter_kernel_sums(scale = scale, nu = 1/2, parent = parent,
                                        daughter = daughter, xlim = xl, ylim = yl)
  d <- sqrt(outer(daughter$x, parent$x, "-")^2 + outer(daughter$y, parent$y, "-")^2)
  expect_equal(sums$sumLogK, sum(log(rowSums(k_exp(d)))), tolerance = 1e-12)
})

test_that("the weights are the ratio of the parent and daughter densities", {
  pats <- patterns()
  theta_0 <- list(kappa = 1, scale = scale, mu = 3)
  theta <- list(kappa = 1.1, scale = 0.33, mu = 2.8)
  log_density <- function(p, th) {
    parent_log_density(pattern = p, kappa = th$kappa) +
      varGamma_daughter_log_density(pattern = p, mu = th$mu, scale = th$scale, nu = 1/2)
  }
  by_hand <- vapply(pats, function(p) exp(log_density(p, theta) - log_density(p, theta_0)), numeric(1))

  w <- importance_sampling_weights_varGamma(kappa = theta$kappa, mu = theta$mu, scale = theta$scale, nu = 1/2,
                                            kappa_0 = theta_0$kappa, mu_0 = theta_0$mu,
                                            scale_0 = theta_0$scale, patternSim = pats)
  expect_equal(w$w_is, by_hand, tolerance = 1e-10)
})

test_that("the cache gives the same weights as computing afresh", {
  pats <- patterns()
  weights <- function(mu, cache = NULL) {
    importance_sampling_weights_varGamma(kappa = 1, mu = mu, scale = 0.33, nu = 1/2,
                                         kappa_0 = 1, mu_0 = 3, scale_0 = scale,
                                         patternSim = pats, daughter_kernel_cache = cache)
  }
  first <- weights(mu = 2.8)
  expect_identical(weights(mu = 3.2, cache = first$daughter_kernel_cache)$w_is, weights(mu = 3.2)$w_is)
})

test_that("the weights are all 1 at the simulation parameters", {
  w <- importance_sampling_weights_varGamma(kappa = 1, mu = 3, scale = scale, nu = 1/2,
                                            kappa_0 = 1, mu_0 = 3, scale_0 = scale,
                                            patternSim = patterns())
  expect_identical(w$w_is, rep(1, 4))
})
