#' Complete estimation procedure using importance sampling + simulation + minimum contrast
#' @description
#' Minimum contrast estimation using trust region approach with simulation based estimation
#' of gradient and SR1 for Hessian approximation.
#' This method also uses importance sampling...
#'
#'
#'
#' @param params Initial guess for parameters to estimate
#' @param parFreeIndex Indices of parameters to estimate
#' @param repRange Repulsion range
#' @param K_hat Estimated K-function for data
#' @param rho_hat Estimated intensity for data
#' @param wq Weights for contrast function
#' @param normalized Normalize IS weights?
#' @param xlims Simulation window x limits
#' @param ylims Simulation window y limits
#' @param nSims Number of simulations
#' @param deltaInit Initial trust radius
#' @param eta Trust region parameter
#' @param deltaMax Greatest trust region radius allowed
#' @param deltaMin Convergence tolerance for trust region radius
#' @param tol Convergence tolerance
#' @param max.iter Maximum iterations
#' @param max.iter.prefit How many iterations to do at most for the pre-fitting
#' @param printProgress Set to TRUE to get updates on progress while running
#'
#' @export
min_contrast_trust_region <- function(params, parFreeIndex, repRange, rho_hat, K_hat,
                                      xlims, ylims, nSims, deltaInit, eta, deltaMax, deltaMin = 0.0001,
                                      wq = c(1000, 1/4), normalized = FALSE, countPrefitLoops = 2, max.iter.prefit = 0,
                                      tol = 10^-8, max.iter = 1000, printProgress = FALSE){
  # TODO: Right now there are multiple while loops where the iteration scheme is ran. These are essentially the same, but with changed parFreeIndex, and wq. Could maybe clean code by making a seperate loop function.
  # Prefit loop: ----
  prefitsCompleted <- 0
  prefitResults <- vector(mode = "list", length = countPrefitLoops)
  while(prefitsCompleted < countPrefitLoops){
    # Pure K-function prefit initialization: ----
    if(printProgress){print("Starting prefitting with only K-function:")}
    tempSteps <- 1
    tempParFreeIndex <- intersect(parFreeIndex, c(1, 2))
    x_sequence_prefit_K <- matrix(data = NA, nrow = max.iter.prefit, ncol = length(tempParFreeIndex))
    x_sequence_prefit_K[1, ] <- params[tempParFreeIndex]
    f_vals_prefit_K <- rep(0, max.iter.prefit)
    update_evaluations_prefit_K <- vector(mode = "list", length = max.iter.prefit)

    simStepRes <- simulation_step(nSims = nSims, params = params, repRange = repRange,
                                  xlims = xlims, ylims = ylims, K_hat = K_hat,
                                  printProgress = printProgress)
    daughter_kernel_cache <- simStepRes$daughter_kernel_cache
    xSim <- params[tempParFreeIndex]
    trust_function <- function(x_new){
      params_new <- params
      params_new[tempParFreeIndex] <- x_new
      res <- contrast_is(simStepRes = simStepRes, params_0 = params,
                         params_new = params_new, rho_hat = rho_hat, K_hat = K_hat,
                         wq = c(0, wq[2]), normalized = normalized, parallel_IS_weights = TRUE,
                         daughter_kernel_cache = daughter_kernel_cache)
      daughter_kernel_cache <<- res$daughter_kernel_cache
      return(res$f_est)
    }

    f_vals_prefit_K[1] <- trust_function(xSim)
    delta <- deltaInit
    converged <- FALSE

    # Pure K-function prefit optimization loop: ----
    while((tempSteps < max.iter.prefit) && !converged){
      if(printProgress){
        print(paste("Starting pure K-function prefit iteration number: ", tempSteps))
      }
      res <- trust_step(x_0 = xSim, delta = delta, trust_function = trust_function)
      params_star <- params
      params_star[tempParFreeIndex] <- res$x_star

      simStepRes_star <- simulation_step(nSims = nSims, params = params_star, repRange = repRange,
                                         xlims = xlims, ylims = ylims, K_hat = K_hat,
                                         printProgress = printProgress)
      f_star <- contrast_is(simStepRes = simStepRes_star, params_0 = params_star,
                            params_new = params_star, rho_hat = rho_hat, K_hat = K_hat,
                            wq = c(0, wq[2]), normalized = normalized, parallel_IS_weights = TRUE)$f_est

      update_eval <- evaluate_improvement(f_old = f_vals_prefit_K[tempSteps], f_new = f_star,
                                          x_new = res$x_star, x_old = x_sequence_prefit_K[tempSteps, ],
                                          delta = delta, deltaMax = deltaMax, eta = eta,
                                          predicted_reduction = f_vals_prefit_K[tempSteps] - res$f_pred)
      update_evaluations_prefit_K[[tempSteps]] <- list(res = res, update_eval = update_eval)

      converged_in_trust_step <- (sum(abs(res$x_star - xSim)) < tol)

      x_sequence_prefit_K[tempSteps + 1, ] <- update_eval$x
      params[tempParFreeIndex] <- update_eval$x
      delta <- update_eval$newDelta
      f_vals_prefit_K[tempSteps + 1] <- update_eval$f_val
      if(update_eval$updated){
        xSim <- update_eval$x
        simStepRes <- simStepRes_star
        daughter_kernel_cache <- simStepRes_star$daughter_kernel_cache

        trust_function <- function(x_new){
          params_new <- params
          params_new[tempParFreeIndex] <- x_new
          res <- contrast_is(simStepRes = simStepRes, params_0 = params,
                             params_new = params_new, rho_hat = rho_hat, K_hat = K_hat,
                             wq = c(0, wq[2]), normalized = normalized, parallel_IS_weights = TRUE,
                             daughter_kernel_cache = daughter_kernel_cache)
          daughter_kernel_cache <<- res$daughter_kernel_cache
          return(res$f_est)
        }
      }
      change_after_step <- sum((x_sequence_prefit_K[tempSteps + 1, ] - x_sequence_prefit_K[tempSteps, ])^2)
      converged <- ((change_after_step < tol) && converged_in_trust_step || (delta < deltaMin))
      tempSteps <- tempSteps + 1
    }
    prefitResults_K <- list(x_seq = x_sequence_prefit_K[1:tempSteps, ],
                            f_vals = f_vals_prefit_K[1:tempSteps],
                            update_evals = update_evaluations_prefit_K[1:tempSteps])
    # Pure mu prefit initialization: ----
    if(printProgress){print("Starting prefitting with only mu:")}
    tempSteps <- 1
    tempParFreeIndex <- c(3)
    x_sequence_prefit_mu <- matrix(data = NA, nrow = max.iter.prefit, ncol = length(tempParFreeIndex))
    x_sequence_prefit_mu[1, ] <- params[tempParFreeIndex]
    f_vals_prefit_mu <- rep(0, max.iter.prefit)
    update_evaluations_prefit_mu <- vector(mode = "list", length = max.iter.prefit)

    simStepRes <- simulation_step(nSims = nSims, params = params, repRange = repRange,
                                  xlims = xlims, ylims = ylims, K_hat = K_hat,
                                  printProgress = printProgress)
    daughter_kernel_cache <- simStepRes$daughter_kernel_cache
    xSim <- params[tempParFreeIndex]
    trust_function <- function(x_new){
      params_new <- params
      params_new[tempParFreeIndex] <- x_new
      res <- contrast_is(simStepRes = simStepRes, params_0 = params,
                         params_new = params_new, rho_hat = rho_hat, K_hat = K_hat,
                         wq = wq, normalized = normalized, parallel_IS_weights = TRUE,
                         daughter_kernel_cache = daughter_kernel_cache)
      daughter_kernel_cache <<- res$daughter_kernel_cache
      return(res$f_est)
    }

    f_vals_prefit_mu[1] <- trust_function(xSim)
    delta <- deltaInit
    converged <- FALSE

    # Pure mu prefit optimization loop: ----
    while((tempSteps < max.iter.prefit) && !converged){
      if(printProgress){
        print(paste("Starting pure mu prefitting iteration number: ", tempSteps))
      }
      res <- trust_step(x_0 = xSim, delta = delta, trust_function = trust_function)
      params_star <- params
      params_star[tempParFreeIndex] <- res$x_star

      simStepRes_star <- simulation_step(nSims = nSims, params = params_star, repRange = repRange,
                                         xlims = xlims, ylims = ylims, K_hat = K_hat,
                                         printProgress = printProgress)
      f_star <- contrast_is(simStepRes = simStepRes_star, params_0 = params_star,
                            params_new = params_star, rho_hat = rho_hat, K_hat = K_hat,
                            wq = wq, normalized = normalized, parallel_IS_weights = TRUE)$f_est

      update_eval <- evaluate_improvement(f_old = f_vals_prefit_mu[tempSteps], f_new = f_star,
                                          x_new = res$x_star, x_old = x_sequence_prefit_mu[tempSteps, ],
                                          delta = delta, deltaMax = deltaMax, eta = eta,
                                          predicted_reduction = f_vals_prefit_mu[tempSteps] - res$f_pred)
      update_evaluations_prefit_mu[[tempSteps]] <- list(res = res, update_eval = update_eval)

      converged_in_trust_step <- (sum(abs(res$x_star - xSim)) < tol)

      x_sequence_prefit_mu[tempSteps + 1, ] <- update_eval$x
      params[tempParFreeIndex] <- update_eval$x
      delta <- update_eval$newDelta
      f_vals_prefit_mu[tempSteps + 1] <- update_eval$f_val
      if(update_eval$updated){
        xSim <- update_eval$x
        simStepRes <- simStepRes_star
        daughter_kernel_cache <- simStepRes_star$daughter_kernel_cache

        trust_function <- function(x_new){
          params_new <- params
          params_new[tempParFreeIndex] <- x_new
          res <- contrast_is(simStepRes = simStepRes, params_0 = params,
                             params_new = params_new, rho_hat = rho_hat, K_hat = K_hat,
                             wq = wq, normalized = normalized, parallel_IS_weights = TRUE,
                             daughter_kernel_cache = daughter_kernel_cache)
          daughter_kernel_cache <<- res$daughter_kernel_cache
          return(res$f_est)
        }
      }
      change_after_step <- sum((x_sequence_prefit_mu[tempSteps + 1, ] - x_sequence_prefit_mu[tempSteps, ])^2)
      converged <- ((change_after_step < tol) && converged_in_trust_step || (delta < deltaMin))
      tempSteps <- tempSteps + 1
    }
    prefitResults_mu <- list(x_seq = x_sequence_prefit_mu[1:tempSteps, ],
                            f_vals = f_vals_prefit_mu[1:tempSteps],
                            update_evals = update_evaluations_prefit_mu[1:tempSteps])

    prefitResults[[prefitsCompleted + 1]] <- list(prefitResults_mu, prefitResults_K)
    prefitsCompleted <- prefitsCompleted + 1
  }

  # Initialize output: ----
  nSteps <- 1
  x_sequence <- matrix(data = NA, nrow = max.iter, ncol = length(parFreeIndex))
  x_sequence[1, ] <- params[parFreeIndex]
  f_vals <- rep(0, max.iter)

  simStepRes <- simulation_step(nSims = nSims, params = params, repRange = repRange,
                                xlims = xlims, ylims = ylims, K_hat = K_hat,
                                printProgress = printProgress)

  daughter_kernel_cache <- simStepRes$daughter_kernel_cache

  f_vals[1] <- trust_function(xSim)
  delta <- deltaInit
  converged <- FALSE

  xSim <- params[parFreeIndex]
  trust_function <- function(x_new){
    params_new <- params
    params_new[parFreeIndex] <- x_new
    res <- contrast_is(simStepRes = simStepRes, params_0 = params,
                       params_new = params_new, rho_hat = rho_hat, K_hat = K_hat,
                       wq = wq, normalized = normalized, parallel_IS_weights = TRUE,
                       daughter_kernel_cache = daughter_kernel_cache)
    daughter_kernel_cache <<- res$daughter_kernel_cache
    return(res$f_est)
  }

  delta <- deltaInit
  converged <- FALSE

  # Optimization loop: ----
  while((nSteps < max.iter) && !converged){
    if(printProgress){
      print(paste("Starting iteration number: ", nSteps))
    }
    res <- trust_step(x_0 = xSim, delta = delta, trust_function = trust_function)
    params_star <- params
    params_star[parFreeIndex] <- res$x_star

    simStepRes_star <- simulation_step(nSims = nSims, params = params_star, repRange = repRange,
                                  xlims = xlims, ylims = ylims, K_hat = K_hat,
                                  printProgress = printProgress)
    f_star <- contrast_is(simStepRes = simStepRes_star, params_0 = params_star,
                          params_new = params_star, rho_hat = rho_hat, K_hat = K_hat,
                          wq = wq, normalized = normalized, parallel_IS_weights = TRUE)$f_est

    update_eval <- evaluate_improvement(f_old = f_vals[nSteps], f_new = f_star,
                                        x_new = res$x_star, x_old = x_sequence[nSteps, ],
                                        delta = delta, deltaMax = deltaMax, eta = eta,
                                        predicted_reduction = f_vals[nSteps] - res$f_pred)

    converged_in_trust_step <- (sum(abs(res$x_star - xSim)) < tol)

    x_sequence[nSteps + 1, ] <- update_eval$x
    params[parFreeIndex] <- update_eval$x
    delta <- update_eval$newDelta
    f_vals[nSteps + 1] <- update_eval$f_val
    if(update_eval$updated){
      xSim <- update_eval$x
      simStepRes <- simStepRes_star
      daughter_kernel_cache <- simStepRes_star$daughter_kernel_cache

      trust_function <- function(x_new){
        params_new <- params
        params_new[parFreeIndex] <- x_new
        res <- contrast_is(simStepRes = simStepRes, params_0 = params,
                           params_new = params_new, rho_hat = rho_hat, K_hat = K_hat,
                           wq = wq, normalized = normalized, parallel_IS_weights = TRUE,
                           daughter_kernel_cache = daughter_kernel_cache)
        daughter_kernel_cache <<- res$daughter_kernel_cache
        return(res$f_est)
      }
    }
    change_after_step <- sum((x_sequence[nSteps + 1, ] - x_sequence[nSteps, ])^2)
    converged <- ((change_after_step < tol) && converged_in_trust_step || (delta < deltaMin))
    nSteps <- nSteps + 1
  }
  return(list(params = x_sequence[1:nSteps, ],
              f_vals = f_vals[1:nSteps]))
}

#' Simulation step in trust region optimization
#' @description
#' This function performs the simulation step in our optimizer.
#'
#' @param nSims Number of simulations
#' @param params params
#' @export
simulation_step <- function(nSims, params, repRange, xlims, ylims, K_hat, printProgress = FALSE){
  if(printProgress){
    print("Simulating pattern: ")
    patternSim <- mcprogress::pmclapply(X = rep(exp(params[1]), nSims), FUN = rThomas_matern_thinned,
                                        scale = exp(params[2]), mu = exp(params[3]),
                                        repulsionRange = repRange, xlims = xlims, ylims = ylims, saveparents = TRUE)
    print("Calculating densities for the simulation parameter:")
    baseline <- baseline_densities(patternSim = patternSim, params = params, printProgress = TRUE)
    print("Estimating baselines: ")
    K_lambda_baseline <- mcprogress::pmclapply(patternSim, estimate_K_lambda_baseline, r_vec = K_hat$r)
  } else {
    patternSim <- parallel::mclapply(X = rep(exp(params[1]), nSims), FUN = rThomas_matern_thinned,
                                     scale = exp(params[2]), mu = exp(params[3]),
                                     repulsionRange = repRange, xlims = xlims, ylims = ylims, saveparents = TRUE)
    K_lambda_baseline <- parallel::mclapply(patternSim, estimate_K_lambda_baseline, r_vec = K_hat$r)
    baseline <- baseline_densities(patternSim = patternSim, params = params, printProgress = FALSE)
  }

  log_f_kappa_0 <- baseline$log_f_kappa_0
  log_fCond_theta_0 <- baseline$log_fCond_theta_0
  daughter_kernel_cache <- baseline$daughter_kernel_cache

  rho_baseline <- sapply(patternSim, estimate_rho_baseline)

  return(list(patternSim = patternSim, rho_baseline = rho_baseline, K_lambda_baseline = K_lambda_baseline,
              log_f_kappa_0 = log_f_kappa_0, log_fCond_theta_0 = log_fCond_theta_0,
              daughter_kernel_cache = daughter_kernel_cache))
}

#' Baseline densities for a freshly simulated ensemble
#' @description
#' Computes the parent and daughter log-densities of `patternSim` under the very
#' parameters it was simulated from, which serve as the denominator of every
#' importance sampling weight taken against that ensemble.
#'
#' The daughter density needs the \eqn{O(n_{daughter} \times n_{parent})} kernel sums at
#' `omega_0`, which are exactly the sums the subsequent trust-region evaluations
#' need whenever they leave omega unchanged. They are therefore returned alongside
#' the densities, ready to seed `daughter_kernel_cache`, so the ensemble is passed
#' over once instead of once here and again on the first evaluation.
#'
#' The parent log-density is `O(1)` per pattern given the enlarged-window area and
#' the parent count, so it is evaluated as a single vectorised expression rather
#' than dispatched per pattern.
#'
#' @param patternSim List of simulated point patterns.
#' @param params Log-scale parameter vector `(log kappa, log omega, log mu)` that
#' `patternSim` was simulated from.
#' @param printProgress Set to TRUE to report progress while computing.
#'
#' @return A list with `log_f_kappa_0`, `log_fCond_theta_0`, and
#' `daughter_kernel_cache` (in the form `importance_sampling_weigths` expects).
baseline_densities <- function(patternSim, params, printProgress = FALSE){
  kappa_0 <- exp(params[1])
  omega_0 <- exp(params[2])
  mu_0    <- exp(params[3])

  B_area   <- vapply(patternSim, function(p) p$B_area, numeric(1))
  n_parent <- vapply(patternSim, function(p) nrow(p$parent), numeric(1))
  log_f_kappa_0 <- B_area*(1 - kappa_0) + n_parent*log(kappa_0)

  if(printProgress){
    kernel_sums_0 <- mcprogress::pmclapply(X = patternSim, FUN = thomas_daughter_kernel_sums,
                                           omega = omega_0)
  } else {
    kernel_sums_0 <- parallel::mclapply(X = patternSim, FUN = thomas_daughter_kernel_sums,
                                        omega = omega_0)
  }
  log_fCond_theta_0 <- vapply(kernel_sums_0, thomas_daughter_log_density_from_sums,
                              numeric(1), mu = mu_0)

  return(list(log_f_kappa_0 = log_f_kappa_0,
              log_fCond_theta_0 = log_fCond_theta_0,
              daughter_kernel_cache = list(omega = omega_0, kernel_sums = kernel_sums_0)))
}

#' Evaluate if there is sufficient improvement
#' @description
#' evaluate improvement
#'
#' @param f_old Value of function at evaluation location
#' @param f_new Value of function at evaluation location
#' @param x_new New parameter values
#' @param x_old Old parameter values
#' @param delta Old trust region radius
#' @param deltaMax Maximal trust region radius
#' @param eta Trust region parameter
#' @param predicted_reduction Reduction in objective with the model function
#'
#' @export
evaluate_improvement <- function(f_old, f_new, x_new, x_old, delta, deltaMax, eta, predicted_reduction){
  if((f_old - f_new)/(predicted_reduction) > 0.75){
    if(sqrt(sum((x_new - x_old)^2)) > 0.99*delta){
      delta <- min(c(2*delta, deltaMax))
    }
  } else if((f_old - f_new)/(predicted_reduction) < 0.1){
    delta <- 0.5*delta
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

  if(!updated && delta > sqrt(sum((x_new - x_old)^2))){
    delta <- min(c(0.8*sqrt(sum((x_new - x_old)^2)), delta))
  }
  return(list(x = x_res, newDelta = delta,
              f_val = f_res, updated = updated))
}


























