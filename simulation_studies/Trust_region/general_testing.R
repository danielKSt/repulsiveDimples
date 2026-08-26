# Load data and libraries: ----
library(spatstat)
library(ggplot2)
library(parallel)
library(mcprogress)
devtools::load_all()

load(file = "simulation_studies/Trust_region/study1.RDa")


# Test the optimizer: ----

# Starting values for parameters from unthinned model:
start_fitted <- thomas.estK(X = K_hat_unions, rmin = 2*min_dist, rmax = 3)
lambda <- -log(1-pi*rho_hat*min_dist^2)/(pi*min_dist^2)
start_fitted <- start_fitted$par

start_params <- c(start_fitted[1], start_fitted[2], lambda/start_fitted[1])
rm(start_fitted, lambda, K_hat_unions, nTrue, sidelength)

res <- min_contrast_trust_region(params = c(log(start_params[1]), log(par_thomas$omega), log(start_params[2])),
                                 par_free_index <- c(1, 3), repRange = min_dist, rho_hat = rho_hat,
                                 K_hat = K_hat, xlims = c(0,3), ylims = c(0,3), nSims = 20000, delta_hat = 0.2,
                                 eta = 0.05, delta_max = 1.0, wq = c(1000, 1/4), normalized = TRUE,
                                 tol = 10^-8, max.iter = 50, printProgress = TRUE)

par_res <- exp(c(res$params[nrow(res$params), 1], log(par_thomas$omega), res$params[nrow(res$params), 2]))
par_thomas
par_res

patternSim <- mcprogress::pmclapply(X = rep(par_res[1], 3000), FUN = rThomas_matern_thinned,
                                    scale = par_res[2], mu = par_res[3],
                                    repulsionRange = min_dist, xlims = c(0,4), ylims = c(0,4), saveparents = TRUE)

rho_baseline <- sapply(X = patternSim, FUN = estimate_rho_baseline)
K_lambda_baseline <- lapply(X = patternSim, FUN = estimate_K_lambda_baseline,
                            r_vec = K_hat$r)

res <- contrast_is(patternSim = patternSim, params_0 = log(par_res),
                   params_new = log(par_res), rho_hat = rho_hat, K_hat = K_hat,
                   rho_baseline = rho_baseline, K_lambda_baseline = K_lambda_baseline,
                   wq = c(1000, 1/4), normalized = TRUE, parallel_IS_weights = TRUE)$f_est

patternSim <- mcprogress::pmclapply(X = rep(par_thomas$kappa, 3000), FUN = rThomas_matern_thinned,
                                    scale = par_res[2], mu =par_thomas$mu,
                                    repulsionRange = min_dist, xlims = c(0,4), ylims = c(0,4), saveparents = TRUE)

rho_baseline <- sapply(X = patternSim, FUN = estimate_rho_baseline)
K_lambda_baseline <- lapply(X = patternSim, FUN = estimate_K_lambda_baseline,
                            r_vec = K_hat$r)

res_true <- contrast_is(patternSim = patternSim, params_0 = log(par_res),
                   params_new = log(par_res), rho_hat = rho_hat, K_hat = K_hat,
                   rho_baseline = rho_baseline, K_lambda_baseline = K_lambda_baseline,
                   wq = c(1000, 1/4), normalized = TRUE, parallel_IS_weights = TRUE)$f_est

trust_function <- function(x_new){
  params_new <- params
  params_new[par_free_index] <- x_new
  res <- contrast_is(patternSim = patternSim, params_0 = params,
                     params_new = params_new, rho_hat = rho_hat, K_hat = K_hat,
                     rho_baseline = rho_baseline, K_lambda_baseline = K_lambda_baseline,
                     wq = wq, normalized = normalized, log_f_kappa_0 = log_f_kappa_0,
                     log_fCond_theta_0 = log_fCond_theta_0, parallel_IS_weights = TRUE,
                     daughter_kernel_cache = daughter_kernel_cache)
  daughter_kernel_cache <<- res$daughter_kernel_cache
  return(res$f_est)
}
