#' Complete estimation procedure using importance sampling + simulation + minimum contrast
#' @description
#' Minimum contrast estimation using a derivative-free trust region approach. Each
#' iteration simulates an ensemble at the current parameters (see
#' \code{\link{simulation_step}}), minimises the contrast function over the trust
#' region by re-weighting that one ensemble with importance sampling (see
#' \code{\link{trust_step}}), and then re-simulates at the candidate to decide whether
#' to accept the step (see \code{\link{evaluate_improvement}}).
#'
#' Before the main fit, the parameters can optionally be pre-fitted in blocks: first
#' `(log kappa, log omega)` against the K-function alone, then `log mu` against the full
#' contrast function. This is repeated `countPrefitLoops` times, and is skipped entirely
#' when `max.iter.prefit` is 0.
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
#' @param max.iter.prefit How many iterations to do at most for the pre-fitting. Set to 0
#' (the default) to skip pre-fitting altogether.
#' @param countPrefitLoops How many times to repeat the pre-fitting cycle, where one cycle
#' is a pure K-function fit of `(log kappa, log omega)` followed by a fit of `log mu`
#' alone. Ignored when `max.iter.prefit` is 0.
#' @param printProgress Set to TRUE to get updates on progress while running
#'
#' @return A list with `params` (the matrix of iterates, one row per step),
#' `f_vals` (the contrast value at each iterate), `update_evals` (the per-step
#' \code{\link{trust_step}} and \code{\link{evaluate_improvement}} output), and
#' `prefitResults` (the same three quantities for each pre-fitting block, or `NULL`
#' when no pre-fitting was done).
#'
#' @export
min_contrast_trust_region <- function(params, parFreeIndex, repRange, rho_hat, K_hat,
                                      xlims, ylims, nSims, deltaInit, eta, deltaMax, deltaMin = 0.0001,
                                      wq = c(1000, 1/4), normalized = FALSE, countPrefitLoops = 2, max.iter.prefit = 10,
                                      tol = 10^-8, max.iter = 1000, printProgress = FALSE){
  if(max.iter < 1){
    stop("max.iter must be at least 1.")
  }

  # Prefit loop: ----
  # One cycle fits (log kappa, log omega) against the K-function alone, then log mu on
  # its own, both restricted to whatever the caller actually left free. Each block is
  # the same trust region iteration as the main fit, only over fewer coordinates and
  # with different contrast weights, so all three are the same call.
  prefitResults <- NULL
  if((max.iter.prefit > 0) && (countPrefitLoops > 0)){
    prefitBlocks <- list(K  = list(index = intersect(parFreeIndex, c(1, 2)), wq = c(0, wq[2])),
                         mu = list(index = intersect(parFreeIndex, 3),       wq = wq))
    prefitResults <- vector(mode = "list", length = countPrefitLoops)

    for(prefitLoop in 1:countPrefitLoops){
      loopResults <- vector(mode = "list", length = length(prefitBlocks))
      names(loopResults) <- names(prefitBlocks)

      for(blockName in names(prefitBlocks)){
        block <- prefitBlocks[[blockName]]
        if(length(block$index) == 0){
          next
        }
        if(printProgress){
          print(paste0("Starting prefit cycle ", prefitLoop, " of ", countPrefitLoops,
                       ", fitting on ", blockName, ":"))
        }
        blockRes <- trust_region_loop(params = params, parFreeIndex = block$index,
                                      repRange = repRange, rho_hat = rho_hat, K_hat = K_hat,
                                      xlims = xlims, ylims = ylims, nSims = nSims,
                                      deltaInit = deltaInit, eta = eta, deltaMax = deltaMax,
                                      deltaMin = deltaMin, wq = block$wq, normalized = normalized,
                                      tol = tol, max.iter = max.iter.prefit,
                                      printProgress = printProgress,
                                      label = paste0("prefit (", blockName, ")"))
        params <- blockRes$params
        loopResults[[blockName]] <- blockRes
      }
      prefitResults[[prefitLoop]] <- loopResults
    }
  }

  # Main fit: ----
  res <- trust_region_loop(params = params, parFreeIndex = parFreeIndex,
                           repRange = repRange, rho_hat = rho_hat, K_hat = K_hat,
                           xlims = xlims, ylims = ylims, nSims = nSims,
                           deltaInit = deltaInit, eta = eta, deltaMax = deltaMax,
                           deltaMin = deltaMin, wq = wq, normalized = normalized,
                           tol = tol, max.iter = max.iter, printProgress = printProgress,
                           label = "main fit")

  return(list(params = res$x_seq,
              f_vals = res$f_vals,
              update_evals = res$update_evals,
              prefitResults = prefitResults))
}

#' One trust region fit over a subset of the parameters
#' @description
#' Runs the trust region iteration to convergence over the coordinates named by
#' `parFreeIndex`, holding the remaining entries of `params` fixed. This is the single
#' loop behind \code{\link{min_contrast_trust_region}}: the pre-fitting blocks and the
#' main fit differ only in which coordinates are free and how the contrast function is
#' weighted, so they are the same call with different `parFreeIndex` and `wq`.
#'
#' Each iteration simulates an ensemble at the current parameters, minimises the
#' importance sampling estimate of the contrast function over the trust region with
#' \code{\link{trust_step}}, re-simulates at the candidate to get an independent
#' estimate there, and accepts or rejects the step with
#' \code{\link{evaluate_improvement}}. Iteration stops once a step both leaves the
#' iterate unchanged and is itself the trust region minimiser, or once the trust radius
#' has shrunk below `deltaMin`.
#'
#' @param params Full log-scale parameter vector `(log kappa, log omega, log mu)` to
#' start from. Entries outside `parFreeIndex` are held fixed throughout.
#' @param parFreeIndex Indices of `params` to estimate.
#' @param repRange Repulsion range.
#' @param rho_hat Estimated intensity for data.
#' @param K_hat Estimated K-function for data.
#' @param xlims Simulation window x limits.
#' @param ylims Simulation window y limits.
#' @param nSims Number of simulations per iteration.
#' @param deltaInit Initial trust radius.
#' @param eta Trust region parameter.
#' @param deltaMax Greatest trust region radius allowed.
#' @param deltaMin Convergence tolerance for trust region radius.
#' @param wq Weights for contrast function.
#' @param normalized Normalize IS weights?
#' @param tol Convergence tolerance.
#' @param max.iter Maximum iterations.
#' @param printProgress Set to TRUE to get updates on progress while running.
#' @param label Name for this fit, used only in the progress messages.
#'
#' @return A list with `params` (the full parameter vector at the final iterate, with the
#' fixed entries untouched), `x_seq` (matrix of free-coordinate iterates, one row per
#' step), `f_vals` (contrast value at each iterate) and `update_evals` (the
#' \code{\link{trust_step}} and \code{\link{evaluate_improvement}} output for each step).
#'
#' @export
trust_region_loop <- function(params, parFreeIndex, repRange, rho_hat, K_hat,
                              xlims, ylims, nSims, deltaInit, eta, deltaMax, deltaMin = 0.0001,
                              wq = c(1000, 1/4), normalized = FALSE,
                              tol = 10^-8, max.iter = 1000, printProgress = FALSE,
                              label = "trust region"){
  # Initialize output: ----
  nSteps <- 1
  x_sequence <- matrix(data = NA, nrow = max.iter, ncol = length(parFreeIndex))
  x_sequence[1, ] <- params[parFreeIndex]
  f_vals <- rep(0, max.iter)
  update_evaluations <- vector(mode = "list", length = max.iter)

  simStepRes <- simulation_step(nSims = nSims, params = params, repRange = repRange,
                                xlims = xlims, ylims = ylims, K_hat = K_hat,
                                printProgress = printProgress)
  daughter_kernel_cache <- simStepRes$daughter_kernel_cache
  xSim <- params[parFreeIndex]

  # `params`, `simStepRes` and `daughter_kernel_cache` are resolved in this frame when
  # the closure runs, not when it is created, so it keeps following them as the loop
  # below reassigns them. It therefore only needs to be defined once.
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

  f_vals[1] <- trust_function(xSim)
  delta <- deltaInit
  converged <- FALSE

  # Optimization loop: ----
  while((nSteps < max.iter) && !converged){
    if(printProgress){
      print(paste0("Starting ", label, " iteration number: ", nSteps))
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
    update_evaluations[[nSteps]] <- list(res = res, update_eval = update_eval)

    converged_in_trust_step <- (sum(abs(res$x_star - xSim)) < tol)

    x_sequence[nSteps + 1, ] <- update_eval$x
    params[parFreeIndex] <- update_eval$x
    delta <- update_eval$newDelta
    f_vals[nSteps + 1] <- update_eval$f_val
    if(update_eval$updated){
      xSim <- update_eval$x
      simStepRes <- simStepRes_star
      daughter_kernel_cache <- simStepRes_star$daughter_kernel_cache
    }
    change_after_step <- sum((x_sequence[nSteps + 1, ] - x_sequence[nSteps, ])^2)
    converged <- (((change_after_step < tol) && converged_in_trust_step) || (delta < deltaMin))
    nSteps <- nSteps + 1
  }

  return(list(params = params,
              x_seq = x_sequence[1:nSteps, , drop = FALSE],
              f_vals = f_vals[1:nSteps],
              update_evals = update_evaluations[seq_len(nSteps - 1)]))
}

#' Simulation step in trust region optimization
#' @description
#' Draws a fresh ensemble of `nSims` thinned Thomas patterns at `params` and computes
#' everything about that ensemble which the subsequent importance sampling evaluations
#' need: the baseline log-densities under `params` itself (the denominator of every
#' importance sampling weight taken against this ensemble), the per-pattern baseline
#' intensity and unnormalized K-function, and the daughter kernel sums at this
#' `params`' omega.
#'
#' Everything in the returned list refers to the same ensemble, so it is passed around
#' as a single object (see \code{\link{contrast_is}}) rather than as separate vectors
#' that have to be kept in step with each other by hand.
#'
#' @param nSims Number of patterns to simulate.
#' @param params Log-scale parameter vector `(log kappa, log omega, log mu)` to simulate
#' from.
#' @param repRange Repulsion range of the Matern II thinning, passed to
#' \code{\link{rThomas_matern_thinned}}.
#' @param xlims Simulation window x limits.
#' @param ylims Simulation window y limits.
#' @param K_hat Estimated K-function for the data. Only `K_hat$r` is used, to fix the
#' radii the baseline K-functions are evaluated on.
#' @param printProgress Set to TRUE to report progress while simulating and computing.
#'
#' @return A list with `patternSim` (the simulated ensemble), `rho_baseline`,
#' `K_lambda_baseline`, `log_f_kappa_0`, `log_fCond_theta_0`, and
#' `daughter_kernel_cache` (in the form \code{\link{importance_sampling_weigths}}
#' expects).
#'
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
