#' One trust region fit over a subset of the parameters
#' @description
#' Runs the trust region iteration to convergence over the coordinates named by
#' `parFreeIndex`, holding the remaining entries of `params` fixed. This is the loop every
#' kernel's optimizer is built on: in \code{\link{min_contrast_trust_region_thomas}} the
#' pre-fitting blocks and the main fit differ only in which coordinates are free and how
#' the contrast function is weighted, so they are the same call with different
#' `parFreeIndex` and `wq`. Everything that depends on the kernel is reached through
#' `kernelType`, in \code{\link{simulation_step}}.
#'
#' Each iteration simulates an ensemble at the current parameters, minimises the
#' importance sampling estimate of the contrast function over the trust region with
#' \code{\link{trust_step}}, re-simulates at the candidate to get an independent
#' estimate there, and accepts or rejects the step with
#' \code{\link{evaluate_improvement}}. Iteration stops once a step both leaves the
#' iterate unchanged and is itself the trust region minimiser, or once the trust radius
#' has shrunk below `deltaMin`.
#'
#' @param params Full log-scale parameter vector to start from, in the order the kernel
#' uses: `(log kappa, log omega, log mu)` for `kernelType = "thomas"`. Entries outside
#' `parFreeIndex` are held fixed throughout.
#' @param parFreeIndex Indices of `params` to estimate.
#' @param kernelType Which cluster process to fit, which decides how the ensembles are
#' simulated and how the importance sampling weights are computed. Currently only
#' `"thomas"`, see \code{\link{simulation_step}}.
#' @param repRange Repulsion range.
#' @param thinningType Type of Matern thinning, 1, 2 or 3, see \code{\link{matern.thinning}}.
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
#'
#' Measured on the parameters rather than on the contrast, the same setting matters much
#' less: 50 paired fits at `nSims = 500` differing only in this argument put the mean L2
#' error in log parameters at 0.625 for 0.9 against 0.663 for 0.99, with 0.99 ahead in 19
#' of the 50 pairs (paired p = 0.15) and costing 1.4 times the wall time. What those runs
#' do settle is where the error comes from. Bias dominates variance about 3:1 in every
#' parameter, and a paired check at `nSims = 5000` put the contrast at 0.0384 at the true
#' parameters against 0.2372 at the median fitted point (15 reps, p = 0.007): the contrast
#' surface prefers the truth, and the optimizer stops short of it. A disappointing fit is
#' therefore evidence about how hard \code{\link{trust_step}} searched -- `max.directions`,
#' and which `method` it used -- rather than about the contrast function or its weights.
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
#' @param carry.directions Keep \code{\link{trust_step}}'s conjugate direction set from
#' one iteration to the next (`TRUE`) instead of rebuilding it from the coordinate axes on
#' every step (`FALSE`, the previous behaviour). Accumulating the set across iterations is
#' the mechanism behind Powell (1964); rebuilding it discards that, and additionally makes
#' every step open with a line search along free parameter 1, which biases the search
#' toward that parameter -- measured on study1's K block, where log kappa overshot its
#' target by 72% while log omega reached only a third of its own. The set is dropped back
#' to the axes whenever it goes singular, since a singular set has lost a search dimension.
#'
#' **Only turn this on together with `max.directions` at `NULL`.** The two interact: with
#' few direction updates per step the carried set refreshes too slowly, a stale conjugate
#' direction persists while the coordinate axes are shuffled out, and the search does worse
#' than if it had restarted. Measured on study1's K block (12 seeds, 10 iterations), median
#' error against the target movement was 0.594 at `max.directions = 1` without carrying,
#' 0.646 at `max.directions = 1` *with* carrying, 0.372 at `NULL` without, and **0.176 at
#' `NULL` with**. Hence the default is `FALSE`, so that the default `max.directions = 1`
#' does not silently land in the worst of those four.
#' @param carry.ess_bounds Keep the parameter values at which the importance sampling was
#' seen to break down from one trust region iteration to the next, when the iteration in
#' between was *rejected*. A rejected step leaves both the iterate and the ensemble
#' untouched, so the region the ensemble supports is exactly the one the previous step
#' mapped, and re-deriving it costs evaluations to learn what is already known. Only the
#' trust radius has changed. The bounds are dropped whenever a step is accepted, since the
#' iterate has moved and a fresh ensemble has been simulated there, and
#' \code{\link{ess_bounds_carried}} additionally refuses any set whose centre does not
#' match the iterate about to be searched from.
#'
#' This is worth much more to `method = "quadratic"` than to `"conjugate"`. The quadratic
#' step opens by bisecting for the breakdown along each coordinate, and carried bounds let
#' those searches discard the part of each bracket beyond a crossing already located
#' without evaluating there; the conjugate step only screens individual line search points,
#' and on a shrunken region most of those were going to be inside the mapped region anyway.
#' Re-running a step at the same centre and ensemble with the radius reduced, over 4
#' ensembles at `nSims = 300`, the evaluations it saved were
#'
#' \tabular{lrr}{
#'   \tab quadratic \tab conjugate \cr
#'   `delta` 0.20 to 0.18 \tab 28-35% \tab 0-2% \cr
#'   `delta` 0.20 to 0.16 \tab 19-38% \tab 1-6% \cr
#'   `delta` 0.20 to 0.10 \tab 0-11%  \tab 0%
#' }
#'
#' The gentler the reduction, the more of the previous step's mapping still applies, and
#' the more there is to save. A step that halves the radius is searching a region small
#' enough that little of what was learned at the old radius bears on it.
#'
#' It costs nothing in accuracy: the iterate was identical with and without the carried
#' bounds in all 24 of those comparisons, and it has to be, since a point the bounds skip
#' is one the search would have evaluated and rejected. What carrying changes is only
#' whether it is paid for.
#' @param method Which \code{\link{trust_step}} method to use, `"quadratic"` (the
#' default) or `"conjugate"` (the original method, and the default until 2026-09-24). See
#' \code{\link{trust_step}}. `subsection_count`, `line_iterations`, `max.directions` and
#' `carry.directions` only apply to the conjugate step.
#' @param interp_fraction How much of the supported box
#' \code{\link{trust_step_quadratic}} spreads its interpolation points over. Ignored when
#' `method` is `"conjugate"`.
#' @param bisection_iterations Bisection steps per ray search inside
#' \code{\link{trust_step_quadratic}}. Ignored when `method` is `"conjugate"`.
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
trust_region_loop <- function(params, parFreeIndex, kernelType = "thomas", repRange, thinningType = 2,
                              rho_hat, K_hat, xlims, ylims, nSims, deltaInit, eta, deltaMax, deltaMin = 0.0001,
                              wq = c(1000, 1/4), normalized = TRUE, eta_trust = 0.6,
                              eta_converged = 0.99, validate.converged = FALSE,
                              subsection_count = 11, line_iterations = 4, max.directions = 1,
                              carry.directions = FALSE, carry.ess_bounds = TRUE,
                              method = c("quadratic", "conjugate"),
                              interp_fraction = 0.25, bisection_iterations = 10,
                              tol = 10^-8, max.iter = 1000, printProgress = FALSE,
                              label = "trust region"){
  method <- match.arg(method)
  # Initialize output: ----
  nSteps <- 1
  x_sequence <- matrix(data = NA, nrow = max.iter, ncol = length(parFreeIndex))
  x_sequence[1, ] <- params[parFreeIndex]
  f_vals <- rep(0, max.iter)
  update_evaluations <- vector(mode = "list", length = max.iter)

  simStepRes <- simulation_step(nSims = nSims, params = params, kernelType = kernelType,
                                repRange = repRange, thinningType = thinningType,
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
  # NULL asks trust_step for the coordinate axes. With carry.directions it is replaced by
  # whatever the step leaves behind, so the conjugate directions accumulate across
  # iterations rather than being rebuilt from the axes every time.
  directions <- NULL
  # Likewise NULL asks trust_step for a fresh set of breakdown bounds. With
  # carry.ess_bounds it is replaced by whatever a rejected step leaves behind, and dropped
  # again the moment a step is accepted and the ensemble moves.
  ess_bounds <- NULL

  # Optimization loop: ----
  while((nSteps < max.iter) && !converged){
    if(printProgress){
      print(paste0("Starting ", label, " iteration number: ", nSteps))
    }
    res <- trust_step(x_0 = xSim, delta = delta, trust_function = trust_function,
                      eta_trust = eta_trust, nSims = nSims, method = method,
                      subsection_count = subsection_count, line_iterations = line_iterations,
                      max.directions = max.directions, directions = directions,
                      interp_fraction = interp_fraction,
                      bisection_iterations = bisection_iterations,
                      ess_bounds = ess_bounds)
    if(carry.directions){
      directions <- res$directions
    }
    if(carry.ess_bounds){
      ess_bounds <- res$ess_bounds
    }
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

    simStepRes_star <- simulation_step(nSims = nSims, params = params_star, kernelType = kernelType,
                                       repRange = repRange, thinningType = thinningType,
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
      # The iterate has moved and the trust function now re-weights a different ensemble,
      # so everything the previous step learned about where that ensemble gives out
      # describes a region this one is no longer searching.
      ess_bounds <- NULL
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
#' Draws a fresh ensemble of `nSims` thinned patterns of the cluster process named by
#' `kernelType` at `params`, and computes everything about that ensemble which the
#' subsequent importance sampling evaluations need: the baseline log-densities under
#' `params` itself (the denominator of every importance sampling weight taken against
#' this ensemble), the per-pattern baseline intensity and unnormalized K-function, and
#' the daughter kernel sums at this `params`' kernel parameters.
#'
#' Everything in the returned list refers to the same ensemble, so it is passed around
#' as a single object (see \code{\link{contrast_is}}) rather than as separate vectors
#' that have to be kept in step with each other by hand. That includes `kernelType`,
#' which \code{\link{contrast_is}} reads to pick the importance sampling weights that
#' match the ensemble.
#'
#' @param nSims Number of patterns to simulate.
#' @param params Log-scale parameter vector to simulate from, `(log kappa, log omega,
#' log mu)` for `kernelType = "thomas"`.
#' @param kernelType Which cluster process to simulate. Currently only `"thomas"`, which
#' simulates with \code{\link{rThomas_matern_thinned}}.
#' @param repRange Repulsion range of the Matern thinning, passed to the simulator.
#' @param thinningType Type of Matern thinning, 1, 2 or 3, passed to the simulator.
#' @param xlims Simulation window x limits.
#' @param ylims Simulation window y limits.
#' @param K_hat Estimated K-function for the data. Only `K_hat$r` is used, to fix the
#' radii the baseline K-functions are evaluated on.
#' @param printProgress Set to TRUE to report progress while simulating and computing.
#'
#' @return A list with `patternSim` (the simulated ensemble), `rho_baseline`,
#' `K_lambda_baseline`, `log_f_kappa_0`, `log_fCond_theta_0`, `daughter_kernel_cache` (in
#' the form the kernel's importance sampling weights expect, see
#' \code{\link{importance_sampling_weights_thomas}}), and `kernelType`.
#'
#' @export
simulation_step <- function(nSims, params, kernelType = "thomas", repRange, thinningType = 2, xlims, ylims,
                            K_hat, printProgress = FALSE){
  simulate_pattern <- switch(kernelType,
                             thomas = function(i){
                               rThomas_matern_thinned(kappa = exp(params[1]), scale = exp(params[2]),
                                                      mu = exp(params[3]), repulsionRange = repRange,
                                                      thinningType = thinningType,
                                                      xlims = xlims, ylims = ylims, saveparents = TRUE)
                             },
                             stop("Unknown kernelType \"", kernelType, "\", use \"thomas\"."))
  if(printProgress){
    print("Simulating pattern: ")
    patternSim <- mcprogress::pmclapply(X = seq_len(nSims), FUN = simulate_pattern)
    print("Calculating densities for the simulation parameter:")
    baseline <- baseline_densities(patternSim = patternSim, params = params, kernelType = kernelType,
                                   printProgress = TRUE)
    print("Estimating baselines: ")
    K_lambda_baseline <- mcprogress::pmclapply(patternSim, estimate_K_lambda_baseline, r_vec = K_hat$r)
  } else {
    patternSim <- parallel::mclapply(X = seq_len(nSims), FUN = simulate_pattern)
    K_lambda_baseline <- parallel::mclapply(patternSim, estimate_K_lambda_baseline, r_vec = K_hat$r)
    baseline <- baseline_densities(patternSim = patternSim, params = params, kernelType = kernelType,
                                   printProgress = FALSE)
  }

  log_f_kappa_0 <- baseline$log_f_kappa_0
  log_fCond_theta_0 <- baseline$log_fCond_theta_0
  daughter_kernel_cache <- baseline$daughter_kernel_cache

  rho_baseline <- sapply(patternSim, estimate_rho_baseline)

  return(list(patternSim = patternSim, rho_baseline = rho_baseline, K_lambda_baseline = K_lambda_baseline,
              log_f_kappa_0 = log_f_kappa_0, log_fCond_theta_0 = log_fCond_theta_0,
              daughter_kernel_cache = daughter_kernel_cache, kernelType = kernelType))
}

#' Baseline densities for a freshly simulated ensemble
#' @description
#' Computes the parent and daughter log-densities of `patternSim` under the very
#' parameters it was simulated from, which serve as the denominator of every
#' importance sampling weight taken against that ensemble. The work is done by the
#' kernel's own function, \code{\link{baseline_densities_thomas}} for `"thomas"`.
#'
#' @param patternSim List of simulated point patterns.
#' @param params Log-scale parameter vector that `patternSim` was simulated from.
#' @param kernelType Which cluster process `patternSim` was simulated from. Currently only
#' `"thomas"`.
#' @param printProgress Set to TRUE to report progress while computing.
#'
#' @return A list with `log_f_kappa_0`, `log_fCond_theta_0`, and `daughter_kernel_cache`.
baseline_densities <- function(patternSim, params, kernelType = "thomas", printProgress = FALSE){
  switch(kernelType,
         thomas = baseline_densities_thomas(patternSim = patternSim, params = params,
                                            printProgress = printProgress),
         stop("Unknown kernelType \"", kernelType, "\", use \"thomas\"."))
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
