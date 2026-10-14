#' Complete estimation procedure for the Variance gamma process using importance sampling + simulation + minimum contrast
#' @description
#' Minimum contrast estimation of a Matern thinned Variance gamma process using a
#' derivative-free trust region approach, see \code{\link{trust_region_loop}}. Each
#' iteration simulates an ensemble at the current parameters (see
#' \code{\link{simulation_step}}), minimises the contrast function over the trust
#' region by re-weighting that one ensemble with importance sampling (see
#' \code{\link{trust_step}}), and then re-simulates at the candidate to decide whether
#' to accept the step (see \code{\link{evaluate_improvement}}).
#'
#' Unlike \code{\link{min_contrast_trust_region_thomas}} there is no pre-fitting: the
#' Thomas pre-fitting rests on the unthinned Thomas K-function not depending on `mu`,
#' which is not expected of the Variance gamma process, so this is a single call to
#' \code{\link{trust_region_loop}} over `parFreeIndex`.
#'
#' @param params Initial guess for the parameters, `(log kappa, log scale, log mu, nu)`.
#' The shape `nu` is on its natural scale, since it can be negative. It is held fixed
#' unless 4 is in `parFreeIndex`.
#' @param parFreeIndex Indices of `params` to estimate. The default leaves `nu` fixed.
#' If `nu` is freed, note that the trust region then mixes a multiplicative restriction
#' on `(kappa, scale, mu)` with an additive one on `nu`, and that
#' \code{\link[spatstat.random]{rVarGamma}} cannot simulate at `nu <= -1/2`.
#' @param thinningType Type of Matern thinning, 1, 2 or 3, see \code{\link{matern.thinning}}.
#' It only changes how the ensembles are simulated. The importance sampling weights are
#' taken over the parents and the unthinned daughters, and the thinning does not depend
#' on `(kappa, scale, mu, nu)`, so the same weights are valid for every type.
#' @param max.directions How many times \code{\link{trust_step}} may update its conjugate
#' direction set. `NULL` (the default, as for the main fit of
#' \code{\link{min_contrast_trust_region_thomas}}) means `2*length(parFreeIndex) - 1`.
#' Only used when `method = "conjugate"`.
#' @inheritParams trust_region_loop
#'
#' @return A list with `params` (the matrix of free-coordinate iterates, one row per
#' step), `f_vals` (the contrast value at each iterate) and `update_evals` (the per-step
#' \code{\link{trust_step}} and \code{\link{evaluate_improvement}} output).
#'
#' @export
min_contrast_trust_region_varGamma <- function(params, parFreeIndex = c(1, 2, 3), repRange, thinningType = 2,
                                               rho_hat, K_hat, xlims, ylims, nSims, deltaInit, eta, deltaMax,
                                               deltaMin = 0.0001, wq = c(1000, 1/4), normalized = TRUE,
                                               eta_trust = 0.6, eta_converged = 0.99,
                                               validate.converged = FALSE,
                                               subsection_count = 11, line_iterations = 4,
                                               max.directions = NULL,
                                               carry.directions = FALSE, carry.ess_bounds = TRUE,
                                               method = c("quadratic", "conjugate"),
                                               interp_fraction = 0.25, bisection_iterations = 10,
                                               tol = 10^-8, max.iter = 1000, printProgress = FALSE){
  if(max.iter < 1){
    stop("max.iter must be at least 1.")
  }
  if(length(params) != 4){
    stop("params must be (log kappa, log scale, log mu, nu), of length 4.")
  }
  method <- match.arg(method)

  res <- trust_region_loop(params = params, parFreeIndex = parFreeIndex,
                           kernelType = "varGamma",
                           repRange = repRange, thinningType = thinningType,
                           rho_hat = rho_hat, K_hat = K_hat,
                           xlims = xlims, ylims = ylims, nSims = nSims,
                           deltaInit = deltaInit, eta = eta, deltaMax = deltaMax,
                           deltaMin = deltaMin, wq = wq, normalized = normalized,
                           eta_trust = eta_trust, eta_converged = eta_converged,
                           validate.converged = validate.converged,
                           subsection_count = subsection_count,
                           line_iterations = line_iterations,
                           max.directions = max.directions,
                           carry.directions = carry.directions,
                           carry.ess_bounds = carry.ess_bounds,
                           method = method,
                           interp_fraction = interp_fraction,
                           bisection_iterations = bisection_iterations,
                           tol = tol, max.iter = max.iter,
                           printProgress = printProgress,
                           label = "main fit")

  return(list(params = res$x_seq,
              f_vals = res$f_vals,
              update_evals = res$update_evals))
}

#' Baseline densities for a freshly simulated Variance gamma ensemble
#' @description
#' Computes the parent and daughter log-densities of `patternSim` under the very
#' parameters it was simulated from, which serve as the denominator of every
#' importance sampling weight taken against that ensemble.
#'
#' The daughter density needs the kernel sums at `(scale_0, nu_0)`, dominated by the
#' window probabilities \eqn{p_W(c)} of every parent (see
#' \code{\link{varGamma_daughter_kernel_sums}}), which are exactly the sums the subsequent
#' trust-region evaluations need whenever they leave `scale` and `nu` unchanged. They are
#' therefore returned alongside the densities, ready to seed `daughter_kernel_cache`, so
#' the ensemble is passed over once instead of once here and again on the first
#' evaluation.
#'
#' @param patternSim List of simulated point patterns.
#' @param params Parameter vector `(log kappa, log scale, log mu, nu)` that `patternSim`
#' was simulated from.
#' @param printProgress Set to TRUE to report progress while computing.
#'
#' @return A list with `log_f_kappa_0`, `log_fCond_theta_0`, and
#' `daughter_kernel_cache` (in the form `importance_sampling_weights_varGamma` expects).
baseline_densities_varGamma <- function(patternSim, params, printProgress = FALSE){
  kappa_0 <- exp(params[1])
  scale_0 <- exp(params[2])
  mu_0    <- exp(params[3])
  nu_0    <- params[4]

  B_area   <- vapply(patternSim, function(p) p$B_area, numeric(1))
  n_parent <- vapply(patternSim, function(p) nrow(p$parent), numeric(1))
  log_f_kappa_0 <- B_area*(1 - kappa_0) + n_parent*log(kappa_0)

  if(printProgress){
    kernel_sums_0 <- mcprogress::pmclapply(X = patternSim, FUN = varGamma_daughter_kernel_sums,
                                           scale = scale_0, nu = nu_0)
  } else {
    kernel_sums_0 <- parallel::mclapply(X = patternSim, FUN = varGamma_daughter_kernel_sums,
                                        scale = scale_0, nu = nu_0)
  }
  log_fCond_theta_0 <- vapply(kernel_sums_0, daughter_log_density_from_sums,
                              numeric(1), mu = mu_0)

  return(list(log_f_kappa_0 = log_f_kappa_0,
              log_fCond_theta_0 = log_fCond_theta_0,
              daughter_kernel_cache = list(scale = scale_0, nu = nu_0, kernel_sums = kernel_sums_0)))
}
