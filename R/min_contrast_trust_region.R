

#' Complete estimation procedure using importance sampling + simulation + minimum contrast
#' @description
#' Minimum contrast estimation using trust region approach with simulation based estimation
#' of gradient and SR1 for Hessian approximation.
#' This method also uses importance sampling...
#'
#'
#' @param tol Convergence tolerance
#' @param max.iter Maximum iterations
#' @param params Initial guess for parameters to estimate
#' @param par_free_index Indices of parameters to estimate
#' @param epsilon Step for derivative calculation
#' @param delta_hat trust radius
#' @param K_hat Estimated K-function for data
#' @param rho_hat Estimated intensity for data
#' @param wqp Weights for contrast function
#' @param normalized Normalize IS weights?
#' @param nSims Number of simulations
#' @param repRange Repulsion range
#' @param xlims Simulation window x limits
#' @param ylims Simulation window y limits
#' @param wqp Weights for contrast function
#' @param eta Trust region parameter
#' @param simulation_threshold simulation_threshold
#' @param r Hyperparameter for SR1 update
#' @param printProgress Set to TRUE to get updates on progress while running
#'
#' @export
minimum_contrast_estimation <- function(params, par_free_index, epsilon, delta_hat,
                                        tol = 10^-8, max.iter = 10000, normalized,
                                        rho_hat, K_hat, nSims, repRange, xlims, ylims, wqp, eta,
                                        simulation_threshold, r, printProgress = FALSE){
  # Initialize output: ----
  nSteps <- 1
  params_sequence <- matrix(data = 0, nrow = max.iter, ncol = length(params))
  f_vals <- rep(0, max.iter)
  gn_vals <- matrix(data = 0, nrow = max.iter, ncol = length(par_free_index))

  # Initaial step ----
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

  evals <- evaluate_contrast_and_gradient(patternSim = patternSim, params = params, params_0 = params,
                                         rho_baseline = rho_baseline, K_lambda_baseline = K_lambda_baseline,
                                         rho_hat = rho_hat, K_hat = K_hat,
                                         par_free_index = par_free_index, epsilon = epsilon,
                                         normalized = normalized, wqp = wqp)

  f_vals[1] <- evals$f_val
  gn_vals[1, ] <- evals$gn_val[par_free_index]

  H <- contrast_hessian_importance_sampling(params = params, epsilon = epsilon,
                                            par_free_index = par_free_index, patternSim = patternSim,
                                            rho_baseline = rho_baseline,
                                            K_lambda_baseline = K_lambda_baseline,
                                            normalized = normalized,
                                            rho_hat = rho_hat, K_hat = K_hat,
                                            wqp = wqp, f_val = f_vals[1])

  trust_region_step <- trust_region_update(evals$gn_val[par_free_index], H, trust_radius = delta_hat, tol = tol)

  params_new <- params
  params_new[par_free_index] <- params[par_free_index] + trust_region_step$params

  update_evaluation <- evaluate_improvement(f_old = f_vals[1], params_new = params_new, params_sim = params_sim,
                                            params_old = params, predicted_reduction = trust_region_step$predicted_red,
                                            patternSim = patternSim, K_hat = K_hat, rho_hat = rho_hat, wqp = wqp,
                                            rho_baseline = rho_baseline, K_lambda_baseline = K_lambda_baseline,
                                            normalized = normalized, delta_hat = delta_hat, eta = eta)

  if(update_evaluation$updateIterate){
    params_sequence[1, ] <- params_new
  } else {
    params_sequence[1, ] <- params
  }

  delta_hat <- update_evaluation$newDelta
  converged <- (gn_vals[1, ]%*%gn_vals[1, ] < tol)

  if(printProgress){
    print(paste("Initial iteration completed using: ", Sys.time() - initTime))
  }

  # Iteration loop ----
  while((nSteps < max.iter) && !converged){
    params <- params_sequence[nSteps, ]
    nSteps <- nSteps + 1
    if(printProgress && (nSteps %% 20 == 0)){
      endTime <- Sys.time()
      print(paste("20 steps completed in: ", endTime - initTime))
      initTime <- Sys.time()
      print(paste("Total steps completed: ", nSteps))
    }
    if(sqrt(sum((params - params_sim)^2)) > simulation_threshold){
      if(printProgress){
        initTime <- Sys.time()
        print(paste("Simulating new pattern at step", nSteps, ": "))
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

      rho_baseline <- sapply(patternSim, estimate_rho_baseline)
      params_sim <- params
    }

    if(update_evaluation$updateIterate){
      evals <- evaluate_contrast_and_gradient(patternSim = patternSim, params = params, params_0 = params_sim,
                                              rho_baseline = rho_baseline, K_lambda_baseline = K_lambda_baseline,
                                              rho_hat = rho_hat, K_hat = K_hat,
                                              par_free_index = par_free_index, epsilon = epsilon,
                                              normalized = normalized, wqp = wqp)

      f_vals[nSteps] <- evals$f_val
      gn_vals[nSteps, ] <- evals$gn_val[par_free_index]
    } else {
      f_vals[nSteps] <- f_vals[nSteps - 1]
      gn_vals[nSteps, ] <- gn_vals[nSteps - 1, ]
    }

    H <- SR1_update(Bk = H, sk = trust_region_step$params, r = r, yk = gn_vals[nSteps, ]- gn_vals[nSteps - 1, ])

    trust_region_step <- trust_region_update(evals$gn_val[par_free_index], H, trust_radius = delta_hat, tol = tol)

    params_new <- params
    params_new[par_free_index] <- params[par_free_index] + trust_region_step$params

    update_evaluation <- evaluate_improvement(f_old = f_vals[nSteps], params_new = params_new, params_sim = params_sim,
                                              params_old = params, predicted_reduction = trust_region_step$predicted_red,
                                              patternSim = patternSim, K_hat = K_hat, rho_hat = rho_hat, wqp = wqp,
                                              rho_baseline = rho_baseline, K_lambda_baseline = K_lambda_baseline,
                                              normalized = normalized, delta_hat = delta_hat, eta = eta)

    if(update_evaluation$updateIterate){
      params_sequence[nSteps, ] <- params_new
    } else {
      params_sequence[nSteps, ] <- params_sequence[nSteps - 1, ]
    }

    delta_hat <- update_evaluation$newDelta
    converged <- (gn_vals[nSteps, ]%*%gn_vals[nSteps, ] < tol)
  }
  return(list(params = params_sequence,
              f_vals = f_vals,
              gradients = gn_vals))
}


#' Evaluate if there is sufficient improvement
#' @description
#' evaluate improvement
#'
#' @param f_old Value of function at evaluation location
#' @param params_new Parameter values at which to approximate the partial derivative
#' @param params_old Parameter values of the previous iteration
#' @param params_sim Parameter values used to simulate  patternSim
#' @param predicted_reduction Step for derivative calculation
#' @param patternSim simulated point pattern to use
#' @param rho_baseline rho_baseline
#' @param K_lambda_baseline K_lambda_baseline
#' @param normalized Normalize IS weights?
#' @param rho_hat Estimated intensity for data
#' @param K_hat Estimated K-function for data
#' @param wqp Weights for contrast function
#' @param delta_hat Old trust region radius
#' @param eta Trust region parameter
#'
#' @export
evaluate_improvement <- function(f_old, params_new, params_old, params_sim, predicted_reduction, patternSim, K_hat, rho_hat, wqp,
                                 rho_baseline, K_lambda_baseline, normalized, delta_hat, eta){
  f_new <- contrast_is(patternSim = patternSim, params_0 = params_sim, params_new = params_new,
                       K_hat = K_hat, rho_hat = rho_hat, wqp = wqp, rho_baseline = rho_baseline,
                       K_lambda_baseline = K_lambda_baseline, normalized = normalized)
  updateIterate <- ((f_old - f_new)/(predicted_reduction) > eta)
  if((f_old - f_new)/(predicted_reduction) > 0.75){
    if(sqrt(sum((params_new - params_old)^2)) > 0.8){
      delta_hat <- 2*delta_hat
    }
  } else if((f_old - f_new)/(predicted_reduction) < 0.1){
    delta_hat <- 0.5*delta_hat
  }
  return(list(updateIterate = updateIterate, newDelta = delta_hat))
}

#' Estimate gradient at some parameter value using importance sampling
#' @description
#' Calculate partial derivative for the k'th parameter
#'
#' @param params Parameter values at which to aproximate the partial derivative
#' @param params_sim Parameter values used in simulating patternSim
#' @param par_free_index Indices of parameters to estimate
#' @param epsilon Step for derivative calculation
#' @param patternSim simulated point pattern to use
#' @param rho_baseline rho_baseline
#' @param K_lambda_baseline K_lambda_baseline
#' @param normalized Normalize IS weights?
#' @param rho_hat Estimated intensity for data
#' @param K_hat Estimated K-function for data
#' @param wqp Weights for contrast function
#' @param f_val Value of function at evaluation location
#'
#' @export
contrast_gradient_importance_sampling <- function(params, params_sim, par_free_index, epsilon, patternSim,
                                                  rho_baseline, K_lambda_baseline, normalized,
                                                  rho_hat, K_hat, wqp, f_val){
  gn <- rep(0, length(params))
  for (i in 1:length(par_free_index)) {
    step_params <- params
    step_params[par_free_index[i]] <- params[par_free_index[i]] + epsilon

    f_fwd <- contrast_is(patternSim = patternSim,
                         params_0 = params_sim,
                         params_new = step_params,
                         rho_hat = rho_hat,
                         K_hat = K_hat,
                         rho_baseline = rho_baseline,
                         K_lambda_baseline = K_lambda_baseline,
                         normalized = normalized,
                         wqp = wqp)

    gn[par_free_index[i]] <- (f_fwd - f_val)/epsilon
  }
  return(gn)
}


#' Estimate Hessian at some parameter value using importance sampling
#' @description
#' See equation 8.21 in Nocedal and Wright.
#'
#' @param params Parameter values at which to aproximate the partial derivative
#' @param par_free_index Indices of parameters to estimate
#' @param epsilon Step for derivative calculation
#' @param patternSim simulated point pattern to use
#' @param rho_baseline rho_baseline
#' @param K_lambda_baseline K_lambda_baseline
#' @param normalized Normalize IS weights?
#' @param rho_hat Estimated intensity for data
#' @param K_hat Estimated K-function for data
#' @param wqp Weights for contrast function
#' @param f_val Value of function at evaluation location
#'
#' @export
contrast_hessian_importance_sampling <- function(params, par_free_index, epsilon, patternSim,
                                                 rho_baseline, K_lambda_baseline, normalized,
                                                 rho_hat, K_hat, wqp, f_val){
  H <- matrix(data = 1, nrow = length(par_free_index), ncol = length(par_free_index))
  for (i in 1:length(par_free_index)) {
    for(j in 1:length(par_free_index)) {
      step_params <- params
      step_params[i] <- params[i] + epsilon
      f_fwd_i <- contrast_is(patternSim = patternSim,
                           params_0 = params,
                           params_new = step_params,
                           rho_hat = rho_hat,
                           K_hat = K_hat,
                           rho_baseline = rho_baseline,
                           K_lambda_baseline = K_lambda_baseline,
                           normalized = normalized,
                           wqp = wqp)

      step_params <- params
      step_params[j] <- params[j] + epsilon
      f_fwd_j <- contrast_is(patternSim = patternSim,
                             params_0 = params,
                             params_new = step_params,
                             rho_hat = rho_hat,
                             K_hat = K_hat,
                             rho_baseline = rho_baseline,
                             K_lambda_baseline = K_lambda_baseline,
                             normalized = normalized,
                             wqp = wqp)

      step_params[i] <- step_params[i] + epsilon
      f_fwd_ij <- contrast_is(patternSim = patternSim,
                             params_0 = params,
                             params_new = step_params,
                             rho_hat = rho_hat,
                             K_hat = K_hat,
                             rho_baseline = rho_baseline,
                             K_lambda_baseline = K_lambda_baseline,
                             normalized = normalized,
                             wqp = wqp)

      H[i,j] <- (f_fwd_ij - f_fwd_i - f_fwd_j + f_val)/epsilon
    }
  }
  return(H)
}



#' Calculate contrast function
#' @description
#' Take the estimated K-function and intensity,
#' together with the K-function and intensity for the parameters,
#' and calculates the contrast fucntion.
#'
#' @param rho_par Intensity for model and parameters
#' @param K_par K-function for model and paramerters
#' @param K_hat Estimated K-function for data
#' @param rho_hat Estimated intensity for data
#' @param wqp Weights for contrast function
#'
#' @export
contrast_function <- function(rho_par, rho_hat, K_par, K_hat,
                              wqp = c(100, 1/4, 2)){

  K_res <- sum(abs(K_par$border^wqp[2] - K_hat$border^wqp[2])^wqp[3])*(K_hat$r[2]-K_hat$r[1])
  rho_res <- abs(rho_par-rho_hat)^2
  return(wqp[1]*rho_res+K_res^(1/wqp[3]))
}


#' Estimate both the contrast function and gradient of the contrast function using importance sampling
#' @description
#' ...?
#'
#' @param patternSim Simulated pattern
#' @param params Parameters to use in estimation
#' @param params_0 Parameters used for simulation
#' @param rho_baseline rho_baseline
#' @param K_lambda_baseline K_lambda_baseline
#' @param rho_hat Estimated intensity for data
#' @param K_hat Estimated K-function for data
#' @param par_free_index Indices of parameters to estimate
#' @param epsilon Step for derivative calculation
#' @param normalized Normalize IS weights?
#' @param wqp Weights for contrast function
#'
#' @export
evaluate_contrast_and_gradient <- function(patternSim, params, params_0, rho_baseline, K_lambda_baseline, rho_hat, K_hat,
                                           par_free_index, epsilon, normalized, wqp){

  if(is.null(params_0)){
    rho_par <- mean(rho_baseline)
    K_par <- K_importance_sampling(w_is = rep(1, length(rho_baseline)), patternSim = NULL,
                                   K_lambda_baseline = K_lambda_baseline,
                                   rho_baseline = rho_baseline)
    f_val <- contrast_function(rho_par = rho_par, rho_hat = rho_hat, K_par = K_par, K_hat = K_hat)
  } else {
    f_val <- contrast_is(patternSim = patternSim, params_0 = params_0, params_new = params,
                         K_hat = K_hat, rho_hat = rho_hat, wqp = wqp, rho_baseline = rho_baseline,
                         K_lambda_baseline = K_lambda_baseline, normalized = normalized)
  }

  gn_val <- contrast_gradient_importance_sampling(params = params, params_sim = params_0, par_free_index = par_free_index,
                                                  epsilon = epsilon, patternSim = patternSim,
                                                  rho_baseline = rho_baseline, K_lambda_baseline = K_lambda_baseline,
                                                  normalized = normalized,
                                                  rho_hat = rho_hat, K_hat = K_hat,
                                                  wqp = wqp, f_val = f_val)

  return(list(f_val = f_val, gn_val = gn_val))
}



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
#' @param wqp Weights for contrast function
#' @param rho_baseline rho_baseline
#' @param K_lambda_baseline K_lambda_baseline
#' @param normalized Normalize IS weights?
#'
#' @export
contrast_is <- function(patternSim, params_0, params_new, rho_hat, K_hat, rho_baseline, K_lambda_baseline, normalized, wqp){
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
                             K_hat = K_hat, K_par = K_est, wqp = wqp)
  return(f_est)
}



















