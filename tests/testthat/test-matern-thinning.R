# The three Matern thinnings differ only in which points of a close pair are removed.
# Run on the same pattern with the same age marks, each type keeps a superset of the
# previous one: a point with no neighbour within range survives type I, a point with
# no older neighbour survives type II, and both of those survive type III.

repulsion <- 0.4

base_pattern <- function(seed = 7) {
  set.seed(seed)
  pp <- spatstat.random::rThomas(kappa = 0.5, scale = 0.3, mu = 10,
                                 win = spatstat.geom::owin(c(0, 10), c(0, 10)))
  data.frame(x = pp$x, y = pp$y)
}

thin <- function(pattern, type, seed = 11) {
  # The same seed gives each type the same age marks.
  set.seed(seed)
  matern.thinning(initialPattern = pattern, repulsionRange = repulsion,
                  xrange = c(0, 10), yrange = c(0, 10), thinningType = type)
}

min_dist <- function(pattern) {
  min(stats::dist(pattern[, c("x", "y")]))
}

keys <- function(pattern) paste(pattern$x, pattern$y)


test_that("every thinning type leaves no pair closer than the repulsion range", {
  pattern <- base_pattern()
  for (type in 1:3) {
    expect_gte(min_dist(thin(pattern, type)), repulsion)
  }
})

test_that("with the same age marks, type I is within type II, which is within type III", {
  pattern <- base_pattern()
  typeI <- keys(thin(pattern, 1))
  typeII <- keys(thin(pattern, 2))
  typeIII <- keys(thin(pattern, 3))

  expect_true(all(typeI %in% typeII))
  expect_true(all(typeII %in% typeIII))
  # On a clustered pattern this dense the three types genuinely differ.
  expect_lt(length(typeI), length(typeII))
  expect_lt(length(typeII), length(typeIII))
})

test_that("type III keeps a point whose only older neighbour was removed", {
  # Oldest to youngest is a, b, c on a line with spacing 0.3: b is removed by a, and
  # c, whose only close neighbour is b, survives type III but not types I or II.
  pattern <- data.frame(x = c(5, 5.3, 5.6), y = c(5, 5, 5))
  set.seed(1)
  ages <- stats::runif(3)
  # matern.thinning draws its own marks, so reorder the points to give a, b, c the
  # ages it is about to draw in decreasing order.
  pattern <- pattern[rank(-ages), ]

  kept <- function(type) sort(thin(pattern, type, seed = 1)$x)
  expect_equal(kept(1), numeric(0))
  expect_equal(kept(2), 5)
  expect_equal(kept(3), c(5, 5.6))
})

test_that("an invalid thinningType is rejected", {
  pattern <- base_pattern()
  expect_error(thin(pattern, 4), "thinningType")
  expect_error(thin(pattern, c(1, 2)), "thinningType")
})

test_that("simulation_step simulates with the thinning type it is given", {
  # Type I is the only type under which no retained point has any unthinned daughter
  # within the repulsion range, so that property shows which type reached the simulator.
  sim <- function(type) {
    old <- options(mc.cores = 1)
    on.exit(options(old), add = TRUE)
    set.seed(3)
    simulation_step(nSims = 4, params = log(c(3, 0.05, 20)), repRange = 0.05,
                    thinningType = type, xlims = c(0, 1), ylims = c(0, 1),
                    K_hat = data.frame(r = seq(0, 0.2, by = 0.05)))$patternSim
  }
  isolated <- function(pattern) {
    if (nrow(pattern$thinned) == 0) return(TRUE)
    d <- spatstat.geom::crossdist(pattern$thinned$x, pattern$thinned$y,
                                  pattern$daughter$x, pattern$daughter$y)
    # Every retained point is itself a daughter, at distance 0.
    all(rowSums(d < 0.05) == 1)
  }

  expect_true(all(vapply(sim(1), isolated, logical(1))))
  expect_false(all(vapply(sim(2), isolated, logical(1))))
})
