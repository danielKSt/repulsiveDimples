# rThomas_matern_thinned and rVarGamma_matern_thinned share their thinning step, and
# the importance sampling code reads the same elements from the patterns of either.

repulsion <- 0.4

sim <- function(process, mu = 5, saveparents = TRUE, seed = 1) {
  set.seed(seed)
  if (process == "Thomas") {
    rThomas_matern_thinned(kappa = 1, scale = 0.3, mu = mu, repulsionRange = repulsion,
                           xlims = c(0, 10), ylims = c(0, 10), saveparents = saveparents)
  } else {
    rVarGamma_matern_thinned(kappa = 1, scale = 0.3, mu = mu, nu = -1/4,
                             repulsionRange = repulsion, xlims = c(0, 10), ylims = c(0, 10),
                             saveparents = saveparents)
  }
}


test_that("both simulators return the same elements", {
  expect_identical(names(sim("VarGamma")), names(sim("Thomas")))
})

test_that("the thinned points lie in the window, no pair closer than the repulsion range", {
  for (process in c("Thomas", "VarGamma")) {
    thinned <- sim(process)$thinned
    expect_gte(min(stats::dist(thinned)), repulsion)
    expect_true(all(thinned$x > 0 & thinned$x < 10 & thinned$y > 0 & thinned$y < 10))
  }
})

test_that("saveparents only changes what is returned, not what is simulated", {
  for (process in c("Thomas", "VarGamma")) {
    expect_identical(sim(process, saveparents = FALSE), sim(process)$thinned)
  }
})

test_that("parents without daughters are kept", {
  # With mu = 0.05 about 95% of the parents have no daughters. spatstat's default,
  # nonempty = TRUE, would drop those and leave about kappa*(1 - exp(-mu)) = 0.05 parents
  # per unit area rather than kappa = 1, and the importance sampling weights assume the
  # parents are a Poisson process of intensity kappa.
  for (process in c("Thomas", "VarGamma")) {
    pattern <- sim(process, mu = 0.05)
    expect_gt(nrow(pattern$parent) / pattern$B_area, 0.7)
  }
})
