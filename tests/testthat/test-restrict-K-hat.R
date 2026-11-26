test_that("r0 = 0 keeps all of K_hat and a positive r0 drops the small radii", {
  r <- seq(from = 0, to = 3, by = 0.05)
  K_hat <- data.frame(r = r, trans = pi * r^2)

  expect_equal(restrict_K_hat(K_hat, 0), K_hat)

  cut <- restrict_K_hat(K_hat, 0.5)
  expect_equal(min(cut$r), 0.5)
  expect_equal(cut$trans, pi * cut$r^2)
})

test_that("restrict_K_hat accepts the fv object K_est.unions returns", {
  set.seed(1)
  pat <- data.frame(x = runif(50, 0, 10), y = runif(50, 0, 10))
  K_hat <- K_est.unions(points_input = list(pat), l = 10, spacing = 5,
                        r_vec = seq(from = 0, to = 2, by = 0.1), timescale = 1)
  cut <- restrict_K_hat(K_hat, 1)
  expect_equal(cut$trans, K_hat$trans[K_hat$r >= 1])
})

test_that("an r0 that leaves fewer than two radii is an error", {
  r <- seq(from = 0, to = 3, by = 0.05)
  expect_error(restrict_K_hat(data.frame(r = r, trans = r), 3), "r0")
})

test_that("a K_hat saved with the old column name border is still accepted", {
  r <- seq(from = 0, to = 3, by = 0.05)
  expect_equal(restrict_K_hat(data.frame(r = r, border = pi * r^2), 0),
               data.frame(r = r, trans = pi * r^2))
})
