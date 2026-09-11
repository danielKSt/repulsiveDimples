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
#' @param eta_trust Smallest effective sample size, as a fraction of `nSims`, an iterate
#' proposed by \code{\link{trust_step}} is allowed to have. Between 0 and 1.
#' @param eta_converged Effective sample size, as a fraction of `nSims`, at or above
#' which the fit is taken to have converged: once the iterate \code{\link{trust_step}}
#' proposes is still this well supported by the ensemble simulated at the current
#' parameters, the importance sampling has not been stretched and there is nothing left
#' to re-simulate for. Expected to be larger than `eta_trust`, and in practice close to 1.
#'
#' A high effective sample size at the proposed iterate says the importance sampling is
#' still reliable there; it does not say the contrast has been minimised, and on study1
#' those come apart. Sweeping this over 0.8/0.9/0.95/0.99/0.999 and scoring each fit's
#' endpoint against a *fresh* 400-pattern ensemble gives mean contrasts of 2.49, 1.85,
#' 1.27, 1.02 and 0.64 against 5.32 at the starting point: the more readily this fires,
#' the worse the fit. Judge it that way and not on the returned `f_vals`, whose last entry
#' at an ESS-converged step is the minimum of the importance sampling surface taken over
#' the same ensemble that defines it, and so is optimistic by up to an order of magnitude
#' (0.143 reported against 2.49 measured, at 0.8).
#' @param validate.converged What to do with the iterate that triggers `eta_converged`.
#' `FALSE` (the default) accepts it straight from the importance sampling surface and
#' stops, never simulating there. `TRUE` puts it through the same fresh simulation and
#' \code{\link{evaluate_improvement}} test as any other step, and only then stops, so the
#' step can still be rejected.
#'
#' The two are worth comparing because they trade different errors. Not validating avoids
#' feeding a second ensemble's sampling error into the decision, but leaves `f_vals`
#' ending on `f_pred`, which is the smallest of several hundred evaluations taken over the
#' ensemble that selected it and is optimistic by roughly an order of magnitude at small
#' `nSims` -- and, measured on study1, does not get better as `nSims` grows, because the
#' line search simply searches harder. Validating costs one simulation and one accept or
#' reject, and puts the last entry of `f_vals` on the same footing as the rest of it.
#'
#' The same footing is not, however, an unbiased one. `evaluate_improvement` accepts a
#' step only when the fresh estimate is low enough, so every entry in `f_vals` has passed
#' a filter that favours a lucky ensemble. On study1 at `nSims = 100` the contrast at one
#' fixed parameter vector varied from 0.12 to 6.01 over twelve independent ensembles, and
#' every value the optimizer reported fell at or below the smallest of those twelve.
#' Neither setting of this argument gives a trustworthy final number, and neither is meant
#' to: `f_vals` is a diagnostic series for watching a fit progress, not a measure of how
#' good the fit is. What the argument changes is which steps get taken. To judge a fit on
#' real data, compare the pair correlation function implied by the fitted parameters
#' against the empirical one, rather than reading anything into the contrast values here.
#' A large spread in the contrast across repeat ensembles at the fitted parameters is
#' worth checking separately: it means `nSims` is too small for the accept/reject decision
#' itself to be meaningful.
#' @param subsection_count.main,subsection_count.prefit Points per line search grid inside
#' \code{\link{trust_step}}, at least 5, for the main fit and for the pre-fitting blocks
#' respectively.
#' @param line_iterations.main,line_iterations.prefit How many times each line search
#' refines its grid, for the main fit and for the pre-fitting blocks. Each round brackets
#' the current minimum two points out on either side, so the searched span shrinks by
#' `4/(subsection_count - 1)` per round and the cost of a line search is
#' `(line_iterations + 1) * subsection_count` evaluations.
#' @param max.directions.main,max.directions.prefit How many times \code{\link{trust_step}}
#' may update its conjugate direction set, for the main fit and for the pre-fitting blocks.
#' `NULL` means the original `2*length(parFreeIndex) - 1` (one fewer than `2*length`,
#' because the original loop counted from 1 and stopped *before* `2*n`).
#'
#' These are split because the two phases want opposite things. Pre-fitting cycles between
#' a `(log kappa, log omega)` block and a `log mu` block, and cycling quickly beats solving
#' either one exactly, so a small value is right there; it also costs nothing on the `mu`
#' block, which has one free parameter and so no direction set to update. The main fit is
#' where the accuracy has to come from: on study1, scoring each fit's endpoint against a
#' fresh 400-pattern ensemble, `max.directions = 1` averaged a contrast of 1.83 against
#' 0.91 for `NULL` -- about twice the error for about a twelfth of the runtime.
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
                                      wq = c(1000, 1/4), normalized = FALSE, eta_trust = 0.6,
                                      eta_converged = 0.99, validate.converged = FALSE,
                                      subsection_count.main = 11, subsection_count.prefit = 11,
                                      line_iterations.main = 4, line_iterations.prefit = 4,
                                      max.directions.main = NULL, max.directions.prefit = 1,
                                      countPrefitLoops = 2, max.iter.prefit = 10,
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
                                      eta_trust = eta_trust, eta_converged = eta_converged,
                                      validate.converged = validate.converged,
                                      subsection_count = subsection_count.prefit,
                                      line_iterations = line_iterations.prefit,
                                      max.directions = max.directions.prefit,
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
                           eta_trust = eta_trust, eta_converged = eta_converged,
                           validate.converged = validate.converged,
                           subsection_count = subsection_count.main,
                           line_iterations = line_iterations.main,
                           max.directions = max.directions.main,
                           tol = tol, max.iter = max.iter,
                           printProgress = printProgress,
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
#' @param eta_trust Smallest effective sample size, as a fraction of `nSims`, an iterate
#' proposed by \code{\link{trust_step}} is allowed to have. Iterates the importance
#' sampling cannot support at that level are rejected inside the trust step rather than
#' being proposed at all.
#' @param eta_converged Effective sample size, as a fraction of `nSims`, at or above
#' which the fit is taken to have converged: once the iterate \code{\link{trust_step}}
#' proposes is still this well supported by the ensemble simulated at the current
#' parameters, the importance sampling has not been stretched and there is nothing left
#' to re-simulate for. Expected to be larger than `eta_trust`, and in practice close to 1.
#'
#' A high effective sample size at the proposed iterate says the importance sampling is
#' still reliable there; it does not say the contrast has been minimised, and on study1
#' those come apart. Sweeping this over 0.8/0.9/0.95/0.99/0.999 and scoring each fit's
#' endpoint against a *fresh* 400-pattern ensemble gives mean contrasts of 2.49, 1.85,
#' 1.27, 1.02 and 0.64 against 5.32 at the starting point: the more readily this fires,
#' the worse the fit. Judge it that way and not on the returned `f_vals`, whose last entry
#' at an ESS-converged step is the minimum of the importance sampling surface taken over
#' the same ensemble that defines it, and so is optimistic by up to an order of magnitude
#' (0.143 reported against 2.49 measured, at 0.8).
#' @param validate.converged What to do with the iterate that triggers `eta_converged`.
#' `FALSE` (the default) accepts it straight from the importance sampling surface and
#' stops, never simulating there. `TRUE` puts it through the same fresh simulation and
#' \code{\link{evaluate_improvement}} test as any other step, and only then stops, so the
#' step can still be rejected.
#'
#' The two are worth comparing because they trade different errors. Not validating avoids
#' feeding a second ensemble's sampling error into the decision, but leaves `f_vals`
#' ending on `f_pred`, which is the smallest of several hundred evaluations taken over the
#' ensemble that selected it and is optimistic by roughly an order of magnitude at small
#' `nSims` -- and, measured on study1, does not get better as `nSims` grows, because the
#' line search simply searches harder. Validating costs one simulation and one accept or
#' reject, and puts the last entry of `f_vals` on the same footing as the rest of it.
#'
#' The same footing is not, however, an unbiased one. `evaluate_improvement` accepts a
#' step only when the fresh estimate is low enough, so every entry in `f_vals` has passed
#' a filter that favours a lucky ensemble. On study1 at `nSims = 100` the contrast at one
#' fixed parameter vector varied from 0.12 to 6.01 over twelve independent ensembles, and
#' every value the optimizer reported fell at or below the smallest of those twelve.
#' Neither setting of this argument gives a trustworthy final number, and neither is meant
#' to: `f_vals` is a diagnostic series for watching a fit progress, not a measure of how
#' good the fit is. What the argument changes is which steps get taken. To judge a fit on
#' real data, compare the pair correlation function implied by the fitted parameters
#' against the empirical one, rather than reading anything into the contrast values here.
#' A large spread in the contrast across repeat ensembles at the fitted parameters is
#' worth checking separately: it means `nSims` is too small for the accept/reject decision
#' itself to be meaningful.
#' @param subsection_count Points per line search grid inside \code{\link{trust_step}},
#' at least 5.
#' @param line_iterations How many times each line search refines its grid.
#' @param max.directions How many times \code{\link{trust_step}} may update its
#' conjugate direction set; `NULL` restores the original `2*length(parFreeIndex) - 1`.
#' The default of 1 trades a slightly cruder individual step for a much cheaper one,
#' which over a whole fit is a win rather than a compromise: the loop simply takes more
#' steps. It suits the pre-fitting blocks in particular, where cycling quickly between the
#' K-function block and the `mu` block beats solving either one exactly. It has no effect
#' on a block with a single free parameter, such as the `mu` block, since that is one line
#' search with no direction set to update.
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
                              wq = c(1000, 1/4), normalized = FALSE, eta_trust = 0.6,
                              eta_converged = 0.99, validate.converged = FALSE,
                              subsection_count = 11, line_iterations = 4, max.directions = 1,
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
    return(list(f_est = res$f_est, ess = res$ess))
  }

  f_vals[1] <- trust_function(xSim)$f_est
  delta <- deltaInit
  converged <- FALSE

  # Optimization loop: ----
  while((nSteps < max.iter) && !converged){
    if(printProgress){
      print(paste0("Starting ", label, " iteration number: ", nSteps))
    }
    res <- trust_step(x_0 = xSim, delta = delta, trust_function = trust_function,
                      eta_trust = eta_trust, nSims = nSims,
                      subsection_count = subsection_count, line_iterations = line_iterations,
                      max.directions = max.directions)
    # The minimiser of the trust function is still well supported by the ensemble simulated
    # at the current iterate, so the importance sampling has not been stretched and the fit
    # has nothing left to re-simulate for.
    ess_converged <- (res$ess_star/nSims > eta_converged)
    if(ess_converged && !validate.converged){
      # Accept that minimiser as it stands. Putting it through evaluate_improvement would
      # mean re-simulating to check a step this ensemble already backs, which feeds a
      # second ensemble's sampling error into the decision. The cost is that `f_vals`
      # ends on `f_pred`, an importance sampling estimate taken over the very ensemble
      # that selected the point, and so optimistic -- badly so at small `nSims`. The
      # earlier entries come from a fresh simulation at each accepted iterate.
      nSteps <- nSteps + 1
      x_sequence[nSteps, ] <- res$x_star
      params[parFreeIndex] <- res$x_star
      f_vals[nSteps] <- res$f_pred
      converged <- TRUE
      break
    }
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
    # `ess_converged` is carried down rather than acted on above, since with
    # `validate.converged` the step still goes through the usual re-simulation and
    # accept/reject first; this assignment would otherwise overwrite that verdict.
    converged <- (ess_converged ||
                    ((change_after_step < tol) && converged_in_trust_step) ||
                    (delta < deltaMin))
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
