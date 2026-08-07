

#' Calculate contrast function
#' @description
#' Take the estimated K-function and intensity,
#' together with the K-function and intensity for the parameters,
#' and calculates the contrast fucntion.
#'
#' @param patternSim Simulated pattern
#' @param params_0 Parameters used in simulation of the point pattern
#' @param params_new Target parameters
#' @param K_hat Estimated K-function for data
#' @param rho_hat Estimated intensity for data
#' @param wq Weights for contrast function
#' @param rho_baseline rho_baseline
#' @param K_lambda_baseline K_lambda_baseline
#' @param normalized Normalize IS weights?
#' @param log_f_kappa_0 Pre-computed density for parent process
#' @param log_fCond_theta_0 Pre-computed density for daughter process
#' @param parallel_IS_weights Set to TRUE to use parallel computing for the importance weights
#'
#' @export
contrast_is <- function(patternSim, params_0, params_new, rho_hat, K_hat,
                        rho_baseline, K_lambda_baseline,
                        log_f_kappa_0 = NULL, log_fCond_theta_0 = NULL,
                        normalized = FALSE, wq = c(1000, 1/4),
                        parallel_IS_weights = TRUE){
  w_is <- importance_sampling_weigths(kappa_0 = exp(params_0[1]),
                                      omega_0 = exp(params_0[2]),
                                      mu_0 = exp(params_0[3]),
                                      kappa = exp(params_new[1]),
                                      omega = exp(params_new[2]),
                                      mu = exp(params_new[3]),
                                      patternSim = patternSim,
                                      log_f_kappa_0 = log_f_kappa_0,
                                      log_fCond_theta_0 = log_fCond_theta_0,
                                      parallel = parallel_IS_weights)
  rho_est <- rho_importance_sampling(w_is = w_is, rho_baseline = rho_baseline, normalized = normalized, patternSim = patternSim)
  K_est <- K_importance_sampling(w_is = w_is, K_lambda_baseline = K_lambda_baseline,
                                 rho_baseline = rho_baseline, normalized = normalized, patternSim = patternSim)
  f_est <- contrast_function(rho_hat = rho_hat, rho_par = rho_est,
                             K_hat = K_hat, K_par = K_est, wq = wq)
  return(f_est)
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

