#' Calculate importance sampling weights for a Variance gamma process
#' @description
#' Calculate importance sampling weights for a Variance gamma process, for moving
#' `(kappa, scale, mu, nu)` away from the `(kappa_0, scale_0, mu_0, nu_0)` the ensemble
#' was simulated with. `nu_0` defaults to `nu`, in which case the shape is not reweighted
#' and `nu` has to be the shape the ensemble was simulated with.
#'
#' @param kappa Kappa parameter
#' @param mu mu parameter
#' @param scale scale parameter, \eqn{\sigma} in \code{\link{varGamma_daughter_log_density}}
#' @param nu Shape parameter of the kernel, greater than -1
#' @param kappa_0 Kappa parameter the ensemble was simulated with
#' @param mu_0 mu parameter the ensemble was simulated with
#' @param scale_0 scale parameter the ensemble was simulated with
#' @param nu_0 Shape parameter the ensemble was simulated with, greater than -1. Defaults
#' to `nu`.
#' @param log_f_kappa_0 Parent null density
#' @param log_fCond_theta_0 Daughter null density
#' @param patternSim List of point patterns to use for estimation
#' @param parallel Set to TRUE to use parallel computing
#' @param daughter_kernel_cache Optional cache from a previous call, as returned in this
#' call's output (a list with `scale`, `nu` and per-pattern `kernel_sums`). The expensive
#' part of the daughter log-density, which depends on `scale` and `nu` but not on `mu`, is
#' only recomputed when `scale` or `nu` differs from the cached value. Pass `NULL` (the
#' default) to always recompute.
#'
#' @return A list with `w_is` (the importance sampling weights) and `daughter_kernel_cache`
#' (to be passed back into the next call for reuse when `scale` and `nu` are unchanged). When the
#' target and simulation parameters coincide the weights are all exactly 1 and
#' `daughter_kernel_cache` is returned unchanged.
#'
#' @export
importance_sampling_weights_varGamma <- function(kappa, mu, scale, nu,
                                                 kappa_0 = NULL, mu_0 = NULL, scale_0 = NULL, nu_0 = nu,
                                                 patternSim, log_f_kappa_0 = NULL, log_fCond_theta_0 = NULL,
                                                 parallel = FALSE, daughter_kernel_cache = NULL){

  if(isTRUE(kappa == kappa_0) && isTRUE(scale == scale_0) && isTRUE(mu == mu_0) && isTRUE(nu == nu_0)){
    return(list(w_is = rep(1, length(patternSim)),
                daughter_kernel_cache = daughter_kernel_cache)) # All parameters being equal gives w_is = 1
  }

  # Parent density
  B_area   <- vapply(patternSim, function(p) p$B_area, numeric(1))
  n_parent <- vapply(patternSim, function(p) nrow(p$parent), numeric(1))

  if(is.null(log_f_kappa_0)){
    log_f_kappa_0 <- B_area*(1 - kappa_0) + n_parent*log(kappa_0)
  }
  log_f_kappa <- B_area*(1 - kappa) + n_parent*log(kappa)

  # Daughters density
  # The expensive O(n_daughter x n_parent) part depends on scale and nu only, so it is
  # computed once per (scale, nu) and reused.
  if(is.null(daughter_kernel_cache) || !isTRUE(all.equal(daughter_kernel_cache$scale, scale)) ||
     !isTRUE(all.equal(daughter_kernel_cache$nu, nu))){
    if(parallel){
      kernel_sums <- parallel::mclapply(X = patternSim, FUN = varGamma_daughter_kernel_sums, scale = scale, nu = nu)
    } else {
      kernel_sums <- lapply(X = patternSim, FUN = varGamma_daughter_kernel_sums, scale = scale, nu = nu)
    }
    daughter_kernel_cache <- list(scale = scale, nu = nu, kernel_sums = kernel_sums)
  }

  log_fCond_theta <- vapply(daughter_kernel_cache$kernel_sums,
                            daughter_log_density_from_sums, numeric(1), mu = mu)

  if(is.null(log_fCond_theta_0)){
    if(isTRUE(all.equal(scale_0, scale)) && isTRUE(all.equal(nu_0, nu))){
      # Baseline and target share the kernel, so the kernel sums just computed are the
      # ones the baseline density needs as well, and only the mu-dependent combine differs.
      log_fCond_theta_0 <- vapply(daughter_kernel_cache$kernel_sums,
                                  daughter_log_density_from_sums, numeric(1), mu = mu_0)
    } else {
      if(parallel){
        kernel_sums_0 <- parallel::mclapply(X = patternSim, FUN = varGamma_daughter_kernel_sums,
                                            scale = scale_0, nu = nu_0)
      } else {
        kernel_sums_0 <- lapply(X = patternSim, FUN = varGamma_daughter_kernel_sums,
                                scale = scale_0, nu = nu_0)
      }
      log_fCond_theta_0 <- vapply(kernel_sums_0, daughter_log_density_from_sums,
                                  numeric(1), mu = mu_0)
    }
  }

  # Return final importance sampling weights
  w_is <- exp(log_f_kappa + log_fCond_theta - log_f_kappa_0 - log_fCond_theta_0)
  return(list(w_is = w_is, daughter_kernel_cache = daughter_kernel_cache))
}

#' Log-density of the daughter pattern given parents for a Variance gamma process
#'
#' @description
#' Computes the log-density of a daughter point pattern \eqn{X} given a parent point
#' pattern \eqn{C}, under a Variance gamma cluster process. Conditional on \eqn{C}, the
#' daughters form a Poisson point process on \eqn{W} with intensity
#' \eqn{\mu \sum_{c_j \in C} k(x - c_j)}, where
#' \deqn{k(u) = \frac{(\|u\|/\sigma)^{\nu} K_{\nu}(\|u\|/\sigma)}{\pi 2^{\nu+1} \sigma^2 \Gamma(\nu+1)}}
#' is the Variance gamma kernel of Jalilian, Guan and Waagepetersen (2013), with
#' \eqn{K_\nu} the modified Bessel function of the second kind. At \eqn{\nu = 1/2} this is
#' the exponential kernel \eqn{k(u) = \exp(-\|u\|/\sigma)/(2\pi\sigma^2)}. The log-density
#' with respect to the unit rate Poisson process is
#' \deqn{
#'   \log f(X \mid C) = |W| - \mu \sum_{c_j \in C} p_W(c_j)
#'   + \sum_{x_i \in X} \log\left(\mu \sum_{c_j \in C} k(x_i - c_j)\right)
#' }
#' where \eqn{p_W(c_j) = \int_W k(x - c_j)\,\mathrm{d}x} is the probability that a
#' daughter of \eqn{c_j} falls inside \eqn{W}, see \code{\link{varGamma_daughter_kernel_sums}}
#' for how it is computed.
#'
#' @param mu       Positive scalar. Mean number of daughters per parent.
#' @param scale    Positive scalar. The scale \eqn{\sigma} of the kernel, `scale` in
#' \code{\link{rVarGamma_matern_thinned}}.
#' @param nu       Scalar greater than -1. The shape \eqn{\nu} of the kernel.
#' @param parent   A \code{data.frame} with columns \code{x} and \code{y} giving parent locations,
#' typically restricted to an enlarged window \eqn{B \supset W}.
#' @param daughter A \code{data.frame} with columns \code{x} and \code{y} giving daughter locations,
#' restricted to the observation window \eqn{W}.
#' @param xlim     Numeric vector of length 2. The x-range of the observation window \eqn{W}.
#' @param ylim     Numeric vector of length 2. The y-range of the observation window \eqn{W}.
#' @param pattern  A list containing elements "parent", "daughter", "xlim_unthinned", and "ylim_unthinned".
#'
#' @return A scalar giving \eqn{\log f(X \mid C)}.
#'
#' @export
varGamma_daughter_log_density <- function(pattern = NULL, mu, scale, nu,
                                          parent = NULL, daughter = NULL,
                                          xlim = NULL, ylim = NULL){
  kernel_sums <- varGamma_daughter_kernel_sums(pattern = pattern, scale = scale, nu = nu,
                                               parent = parent, daughter = daughter,
                                               xlim = xlim, ylim = ylim)
  return(daughter_log_density_from_sums(mu = mu, kernel_sums = kernel_sums))
}

#' Kernel sums for the Variance gamma daughter log-density (the scale and nu dependent part)
#'
#' @description
#' Computes the parts of \code{\link{varGamma_daughter_log_density}} that depend on the
#' kernel but not on \eqn{\mu}: \eqn{\sum_j p_W(c_j)} and
#' \eqn{\sum_i \log \sum_j k(x_i-c_j)}. They are combined with \eqn{\mu} by
#' \code{\link{daughter_log_density_from_sums}}.
#'
#' The 2D Variance gamma kernel does not factorise over the two coordinates, so unlike
#' for the Thomas process \eqn{p_W(c)} has no closed form on a rectangle. It is computed
#' from the normal variance mixture the kernel is simulated from: a daughter is displaced
#' from its parent by \eqn{\sqrt{V} Z}, with \eqn{Z} standard normal in 2D and
#' \eqn{V = 2\sigma^2 G}, \eqn{G \sim \mathrm{Gamma}(\nu + 1, 1)}. Given \eqn{G = g} the
#' probability factorises into normal CDF differences, as for the Thomas process, which
#' leaves one integral per parent,
#' \deqn{p_W(c) = \int_0^\infty P_W(c, \sigma \sqrt{2g})\, \gamma_{\nu+1}(g)\,\mathrm{d}g,}
#' where \eqn{P_W(c, s)} is the probability that a normal displacement with standard
#' deviation \eqn{s} from \eqn{c} lands in \eqn{W} and \eqn{\gamma_{\nu+1}} is the
#' Gamma\eqn{(\nu+1, 1)} density. It is evaluated with \code{\link[stats]{integrate}}
#' after substituting \eqn{g = u^{1/(\nu+1)}}, which turns it into
#' \deqn{p_W(c) = \frac{1}{\Gamma(\nu+2)} \int_0^\infty P_W\left(c, \sigma \sqrt{2u^{1/(\nu+1)}}\right) e^{-u^{1/(\nu+1)}}\,\mathrm{d}u.}
#' For \eqn{\nu < 0} the gamma density is infinite at 0, and for a parent very close to
#' an edge of \eqn{W} the probability also changes sharply there, a combination that
#' made \code{integrate} fail with a roundoff error. The substitution removes the
#' singularity.
#'
#' @param scale    Positive scalar. The scale \eqn{\sigma} of the kernel.
#' @param nu       Scalar greater than -1. The shape \eqn{\nu} of the kernel.
#' @param parent   A \code{data.frame} with columns \code{x} and \code{y} giving parent locations.
#' @param daughter A \code{data.frame} with columns \code{x} and \code{y} giving daughter locations.
#' @param xlim     Numeric vector of length 2. The x-range of the observation window \eqn{W}.
#' @param ylim     Numeric vector of length 2. The y-range of the observation window \eqn{W}.
#' @param pattern  A list containing elements "parent", "daughter", "xlim_unthinned", and "ylim_unthinned".
#'
#' @return A list with \code{area_W}, \code{S_pW} (\eqn{\sum_j p_W(c_j)}),
#' \code{sumLogK} (\eqn{\sum_i \log \sum_j k(x_i-c_j)}), and \code{nDaughter}.
#'
#' @importFrom stats pnorm integrate
#' @export
varGamma_daughter_kernel_sums <- function(pattern = NULL, scale, nu,
                                          parent = NULL, daughter = NULL,
                                          xlim = NULL, ylim = NULL){
  if(!(nu > -1)){
    stop("nu must be greater than -1")
  }
  if(!is.null(pattern)){
    parent <- pattern$parent
    daughter <- pattern$daughter
    xlim <- pattern$xlim_unthinned
    ylim <- pattern$ylim_unthinned
  }
  area_W <- (xlim[2] - xlim[1]) * (ylim[2]-ylim[1])

  S_pW <- 0
  if(nrow(parent) > 0){
    p_W <- mapply(function(cx, cy){
      stats::integrate(function(u){
        g <- u^(1/(nu + 1))
        s <- scale*sqrt(2*g)
        (stats::pnorm(xlim[2], mean = cx, sd = s) - stats::pnorm(xlim[1], mean = cx, sd = s)) *
          (stats::pnorm(ylim[2], mean = cy, sd = s) - stats::pnorm(ylim[1], mean = cy, sd = s)) *
          exp(-g)
      }, lower = 0, upper = Inf, rel.tol = 1e-8)$value
    }, parent$x, parent$y)
    S_pW <- sum(p_W)/gamma(nu + 2)
  }

  sumLogK <- 0
  if (nrow(daughter) > 0) {
    dx <- outer(daughter$x, parent$x, "-")
    dy <- outer(daughter$y, parent$y, "-")
    z <- sqrt(dx*dx + dy*dy)/scale
    # The kernel is z^nu K_nu(z) up to a constant factor, which comes out of the row sums
    # as a single term per daughter.
    E <- z^nu * besselK(z, nu)
    sumLogK <- sum(log(rowSums(E))) - nrow(daughter)*log(pi * 2^(nu + 1) * scale^2 * gamma(nu + 1))
  }

  return(list(area_W = area_W, S_pW = S_pW, sumLogK = sumLogK, nDaughter = nrow(daughter)))
}
