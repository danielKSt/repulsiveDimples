

#' Calculate contrast function
#' @description
#' Take the estimated K-function and intensity,
#' together with the K-function and intensity for the parameters,
#' and calculates the contrast fucntion.
#'
#' @param simStepRes Result from simulation step, see the 'simulation_step' function
#' @param params_0 Parameters used in simulation of the point pattern
#' @param params_new Target parameters
#' @param K_hat Estimated K-function for data
#' @param rho_hat Estimated intensity for data
#' @param wq Weights for contrast function
#' @param normalized Normalize IS weights?
#' @param parallel_IS_weights Set to TRUE to use parallel computing for the importance weights
#' @param daughter_kernel_cache Optional cache passed through to \code{\link{importance_sampling_weigths}};
#' see that function for details. Pass the `daughter_kernel_cache` from this call's return value into
#' the next call to reuse cached work when `params_new`'s omega is unchanged. When `NULL`, the cache
#' `simStepRes` already carries from its own simulation step is used instead of starting cold.
#'
#' @return A list with `f_est` (the contrast value), `daughter_kernel_cache` (to be passed back
#' into the next call for reuse), and `ess` (The effective sample size).
#'
#' @export
contrast_is <- function(simStepRes, params_0, params_new, rho_hat, K_hat,
                        normalized = FALSE, wq = c(1000, 1/4),
                        parallel_IS_weights = TRUE, daughter_kernel_cache = NULL){
  # simStepRes was built by simulation_step(), which already paid for the kernel sums at
  # its own omega. With no cache handed in there is nothing better to start from, and the
  # cache carries the omega it belongs to, so it is simply ignored downstream if
  # params_new moved omega away from it.
  if(is.null(daughter_kernel_cache)){
    daughter_kernel_cache <- simStepRes$daughter_kernel_cache
  }
  is_res <- importance_sampling_weigths(kappa_0 = exp(params_0[1]),
                                        omega_0 = exp(params_0[2]),
                                        mu_0 = exp(params_0[3]),
                                        kappa = exp(params_new[1]),
                                        omega = exp(params_new[2]),
                                        mu = exp(params_new[3]),
                                        patternSim = simStepRes$patternSim,
                                        log_f_kappa_0 = simStepRes$log_f_kappa_0,
                                        log_fCond_theta_0 = simStepRes$log_fCond_theta_0,
                                        parallel = parallel_IS_weights,
                                        daughter_kernel_cache = daughter_kernel_cache)
  w_is <- is_res$w_is
  ess <- (sum(w_is)^2)/sum(w_is^2)
  rho_est <- rho_importance_sampling(w_is = w_is, rho_baseline = simStepRes$rho_baseline, normalized = normalized)
  K_est <- K_importance_sampling(w_is = w_is, K_lambda_baseline = simStepRes$K_lambda_baseline,
                                 rho_baseline = simStepRes$rho_baseline, normalized = normalized)
  f_est <- contrast_function(rho_hat = rho_hat, rho_par = rho_est,
                             K_hat = K_hat, K_par = K_est, wq = wq)
  return(list(f_est = f_est, daughter_kernel_cache = is_res$daughter_kernel_cache, ess = ess))
}

#' Calculate contrast function
#' @description
#' Take the estimated K-function and intensity,
#' together with the K-function and intensity for the parameters,
#' and calculates the contrast fucntion.
#'
#' @param rho_par Intensity for model and parameters
#' @param rho_hat Estimated intensity for data
#' @param K_par K-function for model and paramerters
#' @param K_hat Estimated K-function for data
#' @param wq Weights for contrast function
#'
#' @export
contrast_function <- function(rho_par, rho_hat, K_par, K_hat, wq = c(100, 1/4)){
  K_res <- sum(abs(K_par$border^wq[2] - K_hat$border^wq[2])^2)*(K_hat$r[2]-K_hat$r[1])
  rho_res <- abs(rho_par-rho_hat)^2
  return(wq[1]*rho_res+sqrt(K_res))
}

