#' Calculate importance sampling weights for a Thomas process
#' @description
#' Calculate importance sampling weights for a Thomas process
#'
#'
#' @param kappa Kappa parameter
#' @param mu mu parameter
#' @param omega omega parameter
#' @param kappa_0 Kappa parameter
#' @param mu_0 mu parameter
#' @param omega_0 omega parameter
#' @param log_f_kappa_0 Parent null density
#' @param log_fCond_theta_0 Daughter null density
#' @param patternSim List of point patterns to use for estimation
#' @param parallel Set to TRUE to use parallel computing
#' @param daughter_kernel_cache Optional cache from a previous call, as returned in this
#' call's output (a list with `omega` and per-pattern `kernel_sums`). The expensive,
#' omega-only part of the daughter log-density is only recomputed when `omega` differs
#' from the cached value; otherwise the cached kernel sums are reused. Pass `NULL`
#' (the default) to always recompute.
#'
#' @return A list with `w_is` (the importance sampling weights) and `daughter_kernel_cache`
#' (to be passed back into the next call for reuse when `omega` is unchanged). When the
#' target and simulation parameters coincide the weights are all exactly 1 and
#' `daughter_kernel_cache` is returned unchanged.
#'
#' @export
importance_sampling_weights_thomas <- function(kappa, mu, omega, kappa_0 = NULL, mu_0 = NULL, omega_0 = NULL,
                                               patternSim, log_f_kappa_0 = NULL, log_fCond_theta_0 = NULL,
                                               parallel = FALSE, daughter_kernel_cache = NULL){

  if(isTRUE(kappa == kappa_0) && isTRUE(omega == omega_0) && isTRUE(mu == mu_0)){
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

  # The expensive O(n_daughter x n_parent) part of the daughter density depends
  # on omega only, so it is computed once per omega and reused.
  if(is.null(daughter_kernel_cache) || !isTRUE(all.equal(daughter_kernel_cache$omega, omega))){
    if(parallel){
      kernel_sums <- parallel::mclapply(X = patternSim, FUN = thomas_daughter_kernel_sums, omega = omega)
    } else {
      kernel_sums <- lapply(X = patternSim, FUN = thomas_daughter_kernel_sums, omega = omega)
    }
    daughter_kernel_cache <- list(omega = omega, kernel_sums = kernel_sums)
  }

  log_fCond_theta <- vapply(daughter_kernel_cache$kernel_sums,
                            daughter_log_density_from_sums, numeric(1), mu = mu)

  if(is.null(log_fCond_theta_0)){
    if(isTRUE(all.equal(omega_0, omega))){
      # Baseline and target share the dispersal parameter, so the kernel sums just
      # computed are exactly the ones the baseline density needs - only the cheap
      # mu-dependent combine differs. This is the common case when omega is not
      # one of the parameters being perturbed.
      log_fCond_theta_0 <- vapply(daughter_kernel_cache$kernel_sums,
                                  daughter_log_density_from_sums, numeric(1), mu = mu_0)
    } else {
      if(parallel){
        kernel_sums_0 <- parallel::mclapply(X = patternSim, FUN = thomas_daughter_kernel_sums, omega = omega_0)
      } else {
        kernel_sums_0 <- lapply(X = patternSim, FUN = thomas_daughter_kernel_sums, omega = omega_0)
      }
      log_fCond_theta_0 <- vapply(kernel_sums_0, daughter_log_density_from_sums,
                                  numeric(1), mu = mu_0)
    }
  }

  w_is <- exp(log_f_kappa + log_fCond_theta - log_f_kappa_0 - log_fCond_theta_0)
  return(list(w_is = w_is, daughter_kernel_cache = daughter_kernel_cache))
}

#' Log-density of the daughter pattern given parents for a Thomas cluster process
#'
#' @description
#' Computes the log-density of a daughter point pattern \eqn{X} given a parent
#' point pattern \eqn{C}, under a Thomas cluster process with a Gaussian
#' dispersal kernel. Specifically, conditional on \eqn{C}, the daughters form a
#' Poisson point process on \eqn{W} with conditional intensity
#' \eqn{\lambda(x \mid C) = \mu \sum_{c_j \in C} k(x - c_j)}, where
#' \eqn{k(x - c_j) = \frac{1}{2\pi\omega^2}\exp\!\left(-\frac{\|x-c_j\|^2}{2\omega^2}\right)}
#' is the Gaussian dispersal kernel. The log-density with respect to the unit
#' rate Poisson process reference measure is then:
#' \deqn{
#'   \log f(X \mid C) = |W| - \mu \sum_{c_j \in C} p_W(c_j)
#'   + \sum_{x_i \in X} \log\!\left(\mu \sum_{c_j \in C} k(x_i - c_j)\right)
#' }
#' where \eqn{p_W(c_j) = \int_W k(x - c_j)\,\mathrm{d}x} is the probability
#' that a daughter of parent \eqn{c_j} falls inside \eqn{W}, computed exactly
#' via normal CDF differences for the rectangular window.
#'
#' @param mu       Positive scalar. Mean number of daughters per parent.
#' @param omega    Positive scalar. Standard deviation of the Gaussian dispersal kernel.
#' @param parent   A \code{data.frame} with columns \code{x} and \code{y} giving parent locations,
#' typically restricted to an enlarged window \eqn{B \supset W}.
#' @param daughter A \code{data.frame} with columns \code{x} and \code{y} giving daughter locations,
#' restricted to the observation window \eqn{W}.
#' @param xlim     Numeric vector of length 2. The x-range \eqn{[x_{\min},x_{\max}]} of the observation window \eqn{W}.
#' @param ylim     Numeric vector of length 2. The y-range \eqn{[y_{\min},y_{\max}]} of the observation window \eqn{W}.
#' @param pattern  A list containing elements "parent", "daughter", "xlim", and "ylim" as described above.
#'
#' @return A scalar giving \eqn{\log f(X \mid C)}.
#'
#' @importFrom stats pnorm dnorm
#' @export
thomas_daughter_log_density <- function(pattern = NULL, mu, omega,
                                        parent = NULL, daughter = NULL,
                                        xlim = NULL, ylim = NULL) {
  kernel_sums <- thomas_daughter_kernel_sums(pattern = pattern, omega = omega,
                                             parent = parent, daughter = daughter,
                                             xlim = xlim, ylim = ylim)
  return(daughter_log_density_from_sums(mu = mu, kernel_sums = kernel_sums))
}

#' Kernel sums for the Thomas daughter log-density (the omega-dependent part)
#'
#' @description
#' Computes the parts of \code{\link{thomas_daughter_log_density}} that depend on
#' \eqn{\omega} (and the fixed parent/daughter positions) but not on \eqn{\mu}:
#' \eqn{\sum_j p_W(c_j)} and \eqn{\sum_i \log\left(\sum_j k(x_i-c_j)\right)}. These are
#' the expensive, \eqn{O(n_{daughter} \times n_{parent})} parts of the log-density. Split
#' out so they can be computed once and reused across evaluations that only change
#' \eqn{\mu}, via \code{\link{daughter_log_density_from_sums}}.
#'
#' @param omega    Positive scalar. Standard deviation of the Gaussian dispersal kernel.
#' @param parent   A \code{data.frame} with columns \code{x} and \code{y} giving parent locations.
#' @param daughter A \code{data.frame} with columns \code{x} and \code{y} giving daughter locations.
#' @param xlim     Numeric vector of length 2. The x-range of the observation window \eqn{W}.
#' @param ylim     Numeric vector of length 2. The y-range of the observation window \eqn{W}.
#' @param pattern  A list containing elements "parent", "daughter", "xlim_unthinned", and "ylim_unthinned".
#'
#' @return A list with \code{area_W}, \code{S_pW} (\eqn{\sum_j p_W(c_j)}),
#' \code{sumLogK} (\eqn{\sum_i \log\sum_j k(x_i-c_j)}), and \code{nDaughter}.
#'
#' @importFrom stats pnorm dnorm
#' @export
thomas_daughter_kernel_sums <- function(pattern = NULL, omega,
                                        parent = NULL, daughter = NULL,
                                        xlim = NULL, ylim = NULL){
  if(!is.null(pattern)){
    parent <- pattern$parent
    daughter <- pattern$daughter
    xlim <- pattern$xlim_unthinned
    ylim <- pattern$ylim_unthinned
  }
  area_W <- (xlim[2] - xlim[1]) * (ylim[2]-ylim[1])

  S_pW <- 0
  if(nrow(parent) > 0){
    p_W <- (stats::pnorm(xlim[2], mean = parent$x, sd = omega) - stats::pnorm(xlim[1], mean = parent$x, sd = omega)) *
      (stats::pnorm(ylim[2], mean = parent$y, sd = omega) - stats::pnorm(ylim[1], mean = parent$y, sd = omega))
    S_pW <- sum(p_W)
  }

  sumLogK <- 0
  if (nrow(daughter) > 0) {
    dx <- outer(daughter$x, parent$x, "-")
    dy <- outer(daughter$y, parent$y, "-")
    # dnorm(dx, 0, omega)*dnorm(dy, 0, omega) equals
    # exp(-(dx^2 + dy^2)/(2*omega^2)) / (2*pi*omega^2). Evaluating it in that form
    # needs one exp() over the pair matrix instead of two dnorm() calls, and the
    # constant factor comes out of the row sums as a single term per daughter.
    E <- exp(-(dx*dx + dy*dy) / (2*omega^2))
    sumLogK <- sum(log(rowSums(E))) - nrow(daughter)*log(2*pi*omega^2)
  }

  return(list(area_W = area_W, S_pW = S_pW, sumLogK = sumLogK, nDaughter = nrow(daughter)))
}
