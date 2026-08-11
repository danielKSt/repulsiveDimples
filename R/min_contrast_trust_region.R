

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
                                      xlims, ylims, nSims, delta_hat, eta, delta_max,
                                      wq = c(1000, 1/4), normalized = FALSE,
                                      tol = 10^-8, max.iter = 10000, printProgress = FALSE){
  # Initialize output: ----
  nSteps <- 1
  x_sequence <- matrix(data = NA, nrow = max.iter, ncol = length(par_free_index))
  x_sequence[1, ] <- params[par_free_index]
  f_vals <- rep(0, max.iter)

  if(printProgress){
    initTime <- Sys.time()
    print("Starting initial step")

    print("Simulating pattern: ")
    patternSim <- mcprogress::pmclapply(X = rep(exp(params[1]), nSims), FUN = rThomas_matern_thinned,
                            scale = exp(params[2]), mu = exp(params[3]),
                            repulsionRange = repRange, xlims = xlims, ylims = ylims, saveparents = TRUE)
    print("Calculating densities for the simulation parameter:")
    log_f_kappa_0 <- unlist(mcprogress::pmclapply(X = patternSim, FUN = parent_log_density, kappa = exp(params[1])))
    log_fCond_theta_0 <- unlist(mcprogress::pmclapply(X = patternSim, FUN = thomas_daughter_log_density,
                                               mu = exp(params[3]), omega = exp(params[2])))
    print("Estimating baselines: ")
    K_lambda_baseline <- mcprogress::pmclapply(patternSim, estimate_K_lambda_baseline, r_vec = K_hat$r)
  } else {
    patternSim <- parallel::mclapply(X = rep(exp(params[1]), nSims), FUN = rThomas_matern_thinned,
                                     scale = exp(params[2]), mu = exp(params[3]),
                                     repulsionRange = repRange, xlims = xlims, ylims = ylims, saveparents = TRUE)
    K_lambda_baseline <- parallel::mclapply(patternSim, estimate_K_lambda_baseline, r_vec = K_hat$r)
    log_f_kappa_0 <- unlist(parallel::mclapply(X = patternSim, FUN = parent_log_density, kappa = exp(params[1])))
    log_fCond_theta_0 <- unlist(parallel::mclapply(X = patternSim, FUN = thomas_daughter_log_density,
                                               mu = exp(params[3]), omega = exp(params[2])))
  }

  rho_baseline <- sapply(patternSim, estimate_rho_baseline)

  x_sim <- params[par_free_index]
  daughter_kernel_cache <- NULL
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

  f_vals[1] <- trust_function(x_sim)
  converged <- FALSE

  # Iteration loop ----
  while((nSteps < max.iter) && !converged){
    if(printProgress){print(paste("Starting iteration number: ", nSteps))}
    res <- trust_step(x_0 = x_sim, delta_hat = delta_hat, trust_function = trust_function)
    params_star <- params
    params_star[par_free_index] <- res$x_star
    if(printProgress){
      patternSim_star <- mcprogress::pmclapply(X = rep(exp(params_star[1]), nSims), FUN = rThomas_matern_thinned,
                                               scale = exp(params_star[2]), mu = exp(params_star[3]),
                                               repulsionRange = repRange, xlims = xlims, ylims = ylims, saveparents = TRUE)

      print("Calculating densities for the simulation parameter:")
      log_f_kappa_0_star <- unlist(mcprogress::pmclapply(X = patternSim_star, FUN = parent_log_density, kappa = exp(params_star[1])))
      log_fCond_theta_0_star <- unlist(mcprogress::pmclapply(X = patternSim_star, FUN = thomas_daughter_log_density,
                                                 mu = exp(params_star[3]), omega = exp(params_star[2])))

      print("Estimating baselines: ")
      K_lambda_baseline_star <- mcprogress::pmclapply(patternSim_star, estimate_K_lambda_baseline, r_vec = K_hat$r)
    } else {
      patternSim_star <- parallel::mclapply(X = rep(exp(params_star[1]), nSims), FUN = rThomas_matern_thinned,
                                       scale = exp(params_star[2]), mu = exp(params_star[3]),
                                       repulsionRange = repRange, xlims = xlims, ylims = ylims, saveparents = TRUE)
      K_lambda_baseline_star <- parallel::mclapply(patternSim_star, estimate_K_lambda_baseline, r_vec = K_hat$r)
      log_f_kappa_0_star <- unlist(parallel::mclapply(X = patternSim_star, FUN = parent_log_density, kappa = exp(params_star[1])))
      log_fCond_theta_0_star <- unlist(parallel::mclapply(X = patternSim_star, FUN = thomas_daughter_log_density,
                                              mu = exp(params_star[3]), omega = exp(params_star[2])))
    }
    rho_baseline_star <- sapply(patternSim_star, estimate_rho_baseline)

    f_star <- contrast_is(patternSim = patternSim_star, params_0 = params_star,
                          params_new = params_star, rho_hat = rho_hat, K_hat = K_hat,
                          rho_baseline = rho_baseline_star, K_lambda_baseline = K_lambda_baseline_star,
                          wq = wq, normalized = normalized, log_f_kappa_0 = log_f_kappa_0_star,
                          log_fCond_theta_0 = log_fCond_theta_0_star, parallel_IS_weights = TRUE)$f_est

    update_eval <- evaluate_improvement(f_old = f_vals[nSteps], f_new = f_star,
                                        x_new = res$x_star, x_old = x_sequence[nSteps, ],
                                        delta_hat = delta_hat, delta_max = delta_max, eta = eta,
                                        predicted_reduction = f_vals[nSteps] - res$f_pred)

    x_sequence[nSteps + 1, ] <- update_eval$x
    params[par_free_index] <- update_eval$x
    delta_hat <- update_eval$newDelta
    f_vals[nSteps + 1] <- update_eval$f_val
    if(update_eval$updated){
      patternSim <- patternSim_star
      x_sim <- update_eval$x
      log_f_kappa_0 <- log_f_kappa_0_star
      log_fCond_theta_0 <- log_fCond_theta_0_star
      K_lambda_baseline <- K_lambda_baseline_star
      rho_baseline <- rho_baseline_star
      daughter_kernel_cache <- NULL

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
    }
    change_after_step <- sum((x_sequence[nSteps + 1, ] - x_sequence[nSteps, ])^2)
    converged <- (change_after_step < tol) && (sum(abs(res$x_star - x_sim)) < tol)
    nSteps <- nSteps + 1
  }
  return(list(params = x_sequence[1:nSteps, ],
              f_vals = f_vals[1:nSteps]))
}

#' Do one iteration within trust region
#' @description
#' Iterate from params_sim using importance sampling on the realizations patternSim, simulated from params_sim.
#' The iteration is restricted to the ball centered on the parameters used for simulation,
#' with radius equal to the trust region radius.
#'
#' @param x_0 Parameter values used to simulate patternSim
#' @param delta_hat Trust region radius
#' @param trust_function Trust function
#'
#' @export
trust_step <- function(x_0, delta_hat, trust_function){
  n <- length(x_0)
  if(n == 1){
    x_k <- trust_line_steps(p_k = c(1), x_k = x_0, trust_function = trust_function,
                            bisection_count = 13, delta_hat = delta_hat, x_0 = x_0)
    return(list(x_star = x_k,
                f_pred = trust_function(x_k)))
  }
  p <- matrix(data = 0, nrow = n, ncol = n)
  diag(p) <- 1
  k <- 1
  x_k <- trust_line_steps(x_k = x_0, p_k = p[, 1], x_0 = x_0, delta_hat = delta_hat, trust_function = trust_function)
  hit_boundary <- (sqrt(sum(abs(x_k-x_0)^2)) > 0.99*delta_hat)
  rotated <- FALSE
  while(!hit_boundary && k < 2*n){
    z_j <- matrix(data = 0, nrow = n, ncol = n + 1)
    z_j[, 1] <- x_k
    for(j in 1:n){
      z_j[, j+1] <- trust_line_steps(p_k = p[, j], x_k = z_j[, j], x_0 = x_0,
                                     delta_hat = delta_hat, trust_function = trust_function)
    }
    for(j in 1:(n-1)){
      p[, j] <- p[, j+1]
    }
    p[, n] <- z_j[, n + 1] - z_j[, 1]
    if(det(p)==0){
      if(rotated){
        break
      }
      p <- matrix(data = 0, nrow = n, ncol = n)
      diag(p) <- 1
      p[2, 1] <- -1
      p[1, 2] <- 1
      rotated <- TRUE
    }
    x_k <- trust_line_steps(x_k = z_j[, n + 1], p_k = p[, n], x_0 = x_0,
                            delta_hat = delta_hat, trust_function = trust_function)
    k <- k + 1
    hit_boundary <- (sqrt(sum(abs(x_k-x_0)^2)) > 0.99*delta_hat)
  }

  return(list(x_star = x_k,
              f_pred = trust_function(x_k)))
}

#' Line minimization using interpolation with maximum reach
#' @description
#' Line optimization for the trust steps
#' @param p_k Search direction
#' @param x_k Search starting point
#' @param x_0 Center of trust region
#' @param delta_hat Radius of trust region
#' @param bisection_count Number of bisections
#' @param trust_function Trust function
#'
#' @export
trust_line_steps <- function(p_k, x_k, x_0, delta_hat, trust_function, bisection_count = 11){
  alpha_max <- abs((sum(abs(x_k-x_0)) - delta_hat))
  alpha_k <- seq(from = -alpha_max, to = alpha_max, length.out = bisection_count)
  z_k <- matrix(data = NA, nrow = length(x_k), ncol = bisection_count)
  for(i in 1:length(x_k)){
    z_k[i, ] <- x_k[i] + p_k[i] * alpha_k
  }
  f_vals <- apply(X = z_k, MARGIN = 2, FUN = trust_function)
  j_opt <- which.min(f_vals)
  for(j in 1:5){
    if(j_opt < 3){
      j_opt <- 3
    } else if(j_opt > bisection_count - 2){
      j_opt <- bisection_count - 2
    }
    alpha_k <- seq(from = alpha_k[j_opt-2], to = alpha_k[j_opt+2], length.out = bisection_count)
    for(i in 1:length(x_k)){
      z_k[i, ] <- x_k[i] + p_k[i] * alpha_k
    }
    f_vals <- apply(X = z_k, MARGIN = 2, FUN = trust_function)
    j_opt <- which.min(f_vals)
  }
  return(z_k[ ,j_opt])
}

#' Evaluate if there is sufficient improvement
#' @description
#' evaluate improvement
#'
#' @param f_old Value of function at evaluation location
#' @param f_new Value of function at evaluation location
#' @param x_new New parameter values
#' @param x_old Old parameter values
#' @param delta_hat Old trust region radius
#' @param delta_max Maximal trust region radius
#' @param eta Trust region parameter
#' @param predicted_reduction Reduction in objective with the model function
#'
#' @export
evaluate_improvement <- function(f_old, f_new, x_new, x_old, delta_hat, delta_max, eta, predicted_reduction){
  if((f_old - f_new)/(predicted_reduction) > 0.75){
    if(sqrt(sum((x_new - x_old)^2)) > 0.99*delta_hat){
      delta_hat <- min(c(2*delta_hat, delta_max))
    }
  } else if((f_old - f_new)/(predicted_reduction) < 0.1){
    delta_hat <- 0.5*delta_hat
  }
  if(((f_old - f_new)/predicted_reduction) > eta){
    x_res <- x_new
    f_res <- f_new
    updated <- TRUE
  } else {
    x_res <- x_old
    f_res <- f_old
    updated <- FALSE
  }

  return(list(x = x_res, newDelta = delta_hat,
              f_val = f_res, updated = updated))
}


























