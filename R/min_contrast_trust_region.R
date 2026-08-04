

#' Complete estimation procedure using importance sampling + simulation + minimum contrast
#' @description
#' Minimum contrast estimation using trust region approach with simulation based estimation
#' of gradient and SR1 for Hessian approximation.
#' This method also uses importance sampling...
#'
#'
#'
#' @param params Initial guess for parameters to estimate
#' @param par_free_index Indices of parameters to estimate
#' @param repRange Repulsion range
#' @param K_hat Estimated K-function for data
#' @param rho_hat Estimated intensity for data
#' @param wq Weights for contrast function
#' @param normalized Normalize IS weights?
#' @param xlims Simulation window x limits
#' @param ylims Simulation window y limits
#' @param nSims Number of simulations
#' @param delta_hat trust radius
#' @param eta Trust region parameter
#' @param delta_max Greatest trust region radius allowed
#' @param tol Convergence tolerance
#' @param max.iter Maximum iterations
#' @param printProgress Set to TRUE to get updates on progress while running
#'
#' @export
min_contrast_trust_region <- function(params, par_free_index, repRange, rho_hat, K_hat,
                                      wq, normalized, xlims, ylims, nSims, delta_hat, eta, delta_max,
                                      tol = 10^-8, max.iter = 10000, printProgress = FALSE){
  # Initialize output: ----
  nSteps <- 1
  params_sequence <- matrix(data = 0, nrow = max.iter, ncol = length(params))
  f_vals <- rep(0, max.iter)

  if(printProgress){
    initTime <- Sys.time()
    print("Starting initial step")

    print("Simulating pattern: ")
    patternSim <- mcprogress::pmclapply(X = rep(exp(params[1]), nSims), FUN = rThomas_matern_thinned,
                            scale = exp(params[2]), mu = exp(params[3]),
                            repulsionRange = repRange, xlims = xlims, ylims = ylims, saveparents = TRUE)
    print("Estimating baselines: ")
    K_lambda_baseline <- mcprogress::pmclapply(patternSim, estimate_K_lambda_baseline, r_vec = K_hat$r)
  } else {
    patternSim <- parallel::mclapply(X = rep(exp(params[1]), nSims), FUN = rThomas_matern_thinned,
                                     scale = exp(params[2]), mu = exp(params[3]),
                                     repulsionRange = repRange, xlims = xlims, ylims = ylims, saveparents = TRUE)
    K_lambda_baseline <- parallel::mclapply(patternSim, estimate_K_lambda_baseline, r_vec = K_hat$r)
  }

  params_sim <- params
  rho_baseline <- sapply(patternSim, estimate_rho_baseline)
  f_vals[1] <- contrast_is(patternSim = patternSim, params_0 = params_sim,
                           params_new = params_sim, rho_hat = rho_hat, K_hat = K_hat,
                           rho_baseline = rho_baseline, K_lambda_baseline = K_lambda_baseline)
  converged <- FALSE

  # Iteration loop ----
  while((nSteps < max.iter) && !converged){
    res <- trust_step(f_vals[nSteps], params_sim, patternSim, delta_hat)
    params_star <- res$params
    if(printProgress){
      patternSim_star <- mcprogress::pmclapply(X = rep(exp(params_star[1]), nSims), FUN = rThomas_matern_thinned,
                                               scale = exp(params_star[2]), mu = exp(params_star[3]),
                                               repulsionRange = repRange, xlims = xlims, ylims = ylims, saveparents = TRUE)
      print("Estimating baselines: ")
      K_lambda_baseline_star <- mcprogress::pmclapply(patternSim_star, estimate_K_lambda_baseline, r_vec = K_hat$r)
    } else {
      patternSim_star <- parallel::mclapply(X = rep(exp(params_star[1]), nSims), FUN = rThomas_matern_thinned,
                                       scale = exp(params_star[2]), mu = exp(params_star[3]),
                                       repulsionRange = repRange, xlims = xlims, ylims = ylims, saveparents = TRUE)
      K_lambda_baseline_star <- parallel::mclapply(patternSim_star, estimate_K_lambda_baseline, r_vec = K_hat$r)
    }
    rho_baseline_star <- sapply(patternSim_star, estimate_rho_baseline)

    f_star <- contrast_is(patternSim = patternSim_star, params_0 = params_star,
                          params_new = params_star, rho_hat = rho_hat, K_hat = K_hat,
                          rho_baseline = rho_baseline_star, K_lambda_baseline = K_lambda_baseline_star)

    update_eval <- evaluate_improvement(f_old = f_vals[nSteps], f_new = f_star,
                                        params_new = params_star, params_old = params_sequence[nSteps, ],
                                        delta_hat = delta_hat, delta_max = delta_max, eta = eta,
                                        predicted_reduction = res$predicted_reduction)

    params_sequence[nSteps + 1, ] <- update_eval$params
    delta_hat <- update_eval$newDelta
    f_vals[nSteps + 1] <- update_eval$f_val
    nSteps <- nSteps + 1
  }
  return(list(params = params_sequence,
              f_vals = f_vals))
}

#' Do one iteration within trust region
#' @description
#' Do one
#'
#'
#' @export
trust_step <- function(f_old, patternSim, K_lambda_baseline, rho_baseline, params_sim, delta_hat){

  return(c(f_old, patternSim, K_lambda_baseline, rho_baseline, params_sim, delta_hat))
}


#' Evaluate if there is sufficient improvement
#' @description
#' evaluate improvement
#'
#' @param f_old Value of function at evaluation location
#' @param f_new Value of function at evaluation location
#' @param params_new New parameter values
#' @param params_old Old parameter values
#' @param delta_hat Old trust region radius
#' @param delta_max Maximal trust region radius
#' @param eta Trust region parameter
#' @param predicted_reduction Reduction in objective with the model function
#'
#' @export
evaluate_improvement <- function(f_old, f_new, params_new, params_old, delta_hat, delta_max, eta, predicted_reduction){
  if((f_old - f_new)/(predicted_reduction) > 0.75){
    if(sqrt(sum((params_new - params_old)^2)) > 0.99*delta_hat){
      delta_hat <- min(c(2*delta_hat, delta_max))
    }
  } else if((f_old - f_new)/(predicted_reduction) < 0.1){
    delta_hat <- 0.5*delta_hat
  }
  if(((f_old - f_new)/predicted_reduction) > eta){
    params_res <- params_new
    f_res <- f_new
    updated <- TRUE
  } else {
    params_res <- params_old
    f_res <- f_old
    updated <- FALSE
  }

  return(list(params = params_res, newDelta = delta_hat,
              f_val = f_res, updated = updated))
}


























