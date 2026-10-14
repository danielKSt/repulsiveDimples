# The optimizer pieces for kernelType = "varGamma". Like the Thomas loop tests these
# simulate on a tiny window, so they check bookkeeping and that the ensemble, the
# baseline densities and the weights agree with each other, not statistical quality.

vg_fixture <- function() {
  r <- seq(0, 0.3, by = 0.05)
  list(params = c(log(3), log(0.08), log(4), 1/2),
       repRange = 0.02,
       rho_hat = 12,
       K_hat = data.frame(r = r, border = pi * r^2),
       xlims = c(0, 1), ylims = c(0, 1), nSims = 6)
}

vg_simulation_step <- function(seed = 1) {
  f <- vg_fixture()
  old <- options(mc.cores = 1)
  on.exit(options(old), add = TRUE)
  set.seed(seed)
  simulation_step(nSims = f$nSims, params = f$params, kernelType = "varGamma",
                  repRange = f$repRange, xlims = f$xlims, ylims = f$ylims, K_hat = f$K_hat)
}

test_that("simulation_step's varGamma baseline matches the weights computed from scratch", {
  # The baseline densities and cache that simulation_step stores have to be exactly what
  # importance_sampling_weights_varGamma would compute itself, or every weight taken
  # against the ensemble is off by the difference.
  f <- vg_fixture()
  sim <- vg_simulation_step()
  expect_equal(sim$kernelType, "varGamma")
  expect_equal(sim$daughter_kernel_cache$scale, exp(f$params[2]))
  expect_equal(sim$daughter_kernel_cache$nu, f$params[4])

  new <- f$params + c(0.05, 0.1, -0.05, 0.2)
  args <- list(kappa = exp(new[1]), scale = exp(new[2]), mu = exp(new[3]), nu = new[4],
               kappa_0 = exp(f$params[1]), scale_0 = exp(f$params[2]),
               mu_0 = exp(f$params[3]), nu_0 = f$params[4], patternSim = sim$patternSim)
  w_scratch <- do.call(importance_sampling_weights_varGamma, args)$w_is
  w_stored <- do.call(importance_sampling_weights_varGamma,
                      c(args, list(log_f_kappa_0 = sim$log_f_kappa_0,
                                   log_fCond_theta_0 = sim$log_fCond_theta_0)))$w_is
  expect_equal(w_stored, w_scratch, tolerance = 1e-10)
})

test_that("contrast_is picks the varGamma weights from the ensemble", {
  f <- vg_fixture()
  sim <- vg_simulation_step()
  # At the simulation parameters every weight is 1, so the ESS is the ensemble size.
  same <- contrast_is(simStepRes = sim, params_0 = f$params, params_new = f$params,
                      rho_hat = f$rho_hat, K_hat = f$K_hat, normalized = TRUE,
                      parallel_IS_weights = FALSE)
  expect_equal(same$ess, f$nSims)
  # Moving only mu reuses the stored kernel sums, and has to give the same answer as
  # recomputing them.
  new <- f$params + c(0, 0, 0.1, 0)
  cached <- contrast_is(simStepRes = sim, params_0 = f$params, params_new = new,
                        rho_hat = f$rho_hat, K_hat = f$K_hat, normalized = TRUE,
                        parallel_IS_weights = FALSE)
  cold <- contrast_is(simStepRes = sim, params_0 = f$params, params_new = new,
                      rho_hat = f$rho_hat, K_hat = f$K_hat, normalized = TRUE,
                      parallel_IS_weights = FALSE,
                      daughter_kernel_cache = list(scale = -1, nu = -1, kernel_sums = NULL))
  expect_equal(cached$f_est, cold$f_est, tolerance = 1e-10)
  expect_true(cached$ess < f$nSims)
})

test_that("min_contrast_trust_region_varGamma runs and leaves nu fixed by default", {
  f <- vg_fixture()
  old <- options(mc.cores = 1)
  on.exit(options(old), add = TRUE)
  set.seed(4)
  res <- min_contrast_trust_region_varGamma(params = f$params, repRange = f$repRange,
                                            rho_hat = f$rho_hat, K_hat = f$K_hat,
                                            xlims = f$xlims, ylims = f$ylims,
                                            nSims = f$nSims, deltaInit = 0.2, eta = 0.05,
                                            deltaMax = 1, wq = c(1, 1/4), eta_trust = 0.5,
                                            max.iter = 3, bisection_iterations = 4)
  expect_equal(ncol(res$params), 3L)
  expect_equal(nrow(res$params), length(res$f_vals))
  expect_true(all(is.finite(res$f_vals)))
})

test_that("min_contrast_trust_region_varGamma refuses a Thomas-length parameter vector", {
  f <- vg_fixture()
  expect_error(min_contrast_trust_region_varGamma(params = f$params[1:3], repRange = f$repRange,
                                                  rho_hat = f$rho_hat, K_hat = f$K_hat,
                                                  xlims = f$xlims, ylims = f$ylims,
                                                  nSims = f$nSims, deltaInit = 0.2,
                                                  eta = 0.05, deltaMax = 1),
               "length 4")
})
