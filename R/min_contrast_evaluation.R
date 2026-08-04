

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
#'
#' @export
contrast_is <- function(patternSim, params_0, params_new, rho_hat, K_hat,
                        rho_baseline, K_lambda_baseline,
                        normalized = FALSE, wq = c(1000, 1/4)){
  w_is <- importance_sampling_weigths(kappa_0 = exp(params_0[1]),
                                      mu_0 = exp(params_0[2]),
                                      omega_0 = exp(params_0[3]),
                                      kappa = exp(params_new[1]),
                                      mu = exp(params_new[2]),
                                      omega = exp(params_new[3]),
                                      patternSim = patternSim)
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

