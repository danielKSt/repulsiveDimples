#' Use importance sampling to estimate the K-function
#' @description
#' Use importance sampling to estimate the Ripleys' K-function for parameter values
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
#' (to be passed back into the next call for reuse when `omega` is unchanged).
#'
#' @export
importance_sampling_weigths <- function(kappa, mu, omega, kappa_0 = NULL, mu_0 = NULL, omega_0 = NULL,
                                        patternSim, log_f_kappa_0 = NULL, log_fCond_theta_0 = NULL,
                                        parallel = FALSE, daughter_kernel_cache = NULL){
  if(is.null(log_f_kappa_0)){
    if(parallel){
      log_f_kappa_0 <- unlist(parallel::mclapply(X = patternSim, FUN = parent_log_density, kappa = kappa_0))
    } else {
      log_f_kappa_0 <- sapply(X = patternSim, FUN = parent_log_density, kappa = kappa_0)
    }
  }
  if(is.null(log_fCond_theta_0)){
    if(parallel){
      log_fCond_theta_0 <- unlist(parallel::mclapply(X = patternSim, FUN = thomas_daughter_log_density,
                                              mu = mu_0, omega = omega_0))
    } else{
      log_fCond_theta_0 <- sapply(X = patternSim, FUN = thomas_daughter_log_density,
                                  mu = mu_0, omega = omega_0)
    }
  }

  if(is.null(daughter_kernel_cache) || !isTRUE(all.equal(daughter_kernel_cache$omega, omega))){
    if(parallel){
      kernel_sums <- parallel::mclapply(X = patternSim, FUN = thomas_daughter_kernel_sums, omega = omega)
    } else {
      kernel_sums <- lapply(X = patternSim, FUN = thomas_daughter_kernel_sums, omega = omega)
    }
    daughter_kernel_cache <- list(omega = omega, kernel_sums = kernel_sums)
  }

  if(parallel){
    log_f_kappa <- unlist(parallel::mclapply(X = patternSim, FUN = parent_log_density, kappa = kappa))
  } else {
    log_f_kappa <- sapply(X = patternSim, FUN = parent_log_density, kappa = kappa)
  }
  log_fCond_theta <- sapply(X = daughter_kernel_cache$kernel_sums,
                            FUN = thomas_daughter_log_density_from_sums, mu = mu)

  w_is <- exp(log_f_kappa + log_fCond_theta - log_f_kappa_0 - log_fCond_theta_0)
  return(list(w_is = w_is, daughter_kernel_cache = daughter_kernel_cache))
}


#' Use importance sampling to estimate the K-function
#' @description
#' Use importance sampling to estimate the Ripleys' K-function for parameter values
#'
#' @param w_is Importance sampling weights
#' @param patternSim List of point patterns to use for estimation
#' @param K_lambda_baseline K_lambda_baseline
#' @param r_vec r_vec
#' @param rho_baseline rho_baseline
#' @param normalized Set to TRUE if you want to normalize the importance sampling ratios
#'
#' @export
K_importance_sampling <- function(w_is, K_lambda_baseline, rho_baseline,
                                  r_vec = NULL, patternSim = NULL, normalized = FALSE){

  K_lambda <- K_lambda_importance_sampling(w_is, K_lambda_baseline, r_vec, patternSim, normalized)
  lambda_squared <- rho_importance_sampling(w_is, patternSim = patternSim,
                                            rho_baseline = rho_baseline, normalized = normalized)^2

  K_res <- K_lambda
  K_res$border <- K_res$border/lambda_squared
  return(K_res)
}

#' Use importance sampling to estimate the "unnormalized" K-function
#' @description
#' Use importance sampling to estimate the Ripleys' K-function for parameter values
#'
#' @param w_is Importance sampling weights
#' @param patternSim List of point patterns to use for estimation
#' @param K_lambda_baseline K_lambda_baseline
#' @param r_vec r_vec
#' @param normalized Set to TRUE if you want to normalize the importance sampling ratios
#'
#' @export
K_lambda_importance_sampling <- function(w_is, K_lambda_baseline, r_vec = NULL,
                                         patternSim = NULL, normalized = FALSE){
  if(is.null(K_lambda_baseline)){
    K_lambda_baseline <- lapply(patternSim, estimate_K_lambda_baseline, r_vec = r_vec)
  }

  K_res <- K_lambda_baseline[[1]]
  if(normalized){
    w_is <- w_is/sum(w_is)
  }

  K_res$border <- K_lambda_baseline[[1]]$border*w_is[1]
  for (i in 2:length(K_lambda_baseline)) {
    K_res$border <- K_res$border + K_lambda_baseline[[i]]$border*w_is[i]
  }
  if(!normalized){
    K_res$border <- K_res$border/length(K_lambda_baseline)
  }
  return(K_res)
}

#' Estimate the "unnormalized" K-function directly for a given point pattern
#' @description
#' Estimate the "unnormalized" K-function directly for a given point pattern, this function does only the "border" correction.
#' Pairwise distances and border-correction weights are computed once and reused for every radius in `r_vec`,
#' via a sorted cumulative sum, rather than recomputing the full pairwise arrays per radius.
#'
#' @param pattern The point pattern to estimate for
#' @param r_vec Vector of radii to use
#'
#' @export
estimate_K_lambda_baseline <- function(pattern, r_vec){
  x     <- pattern$thinned$x
  y     <- pattern$thinned$y
  n     <- nrow(pattern$thinned)
  xlim1 <- pattern$xlim_thinned[1]
  xlim2 <- pattern$xlim_thinned[2]
  ylim1 <- pattern$ylim_thinned[1]
  ylim2 <- pattern$ylim_thinned[2]

  if(n < 2){
    return(data.frame(r = r_vec, border = 0))
  }

  dx <- outer(x, x, "-")
  dy <- outer(y, y, "-")
  ut <- upper.tri(dx)
  dx <- dx[ut]
  dy <- dy[ut]

  d <- sqrt(dx^2 + dy^2)
  w <- 1 / ((xlim2 - xlim1 - abs(dx)) * (ylim2 - ylim1 - abs(dy)))

  ord <- order(d)
  d <- d[ord]
  cw <- cumsum(w[ord])

  idx <- findInterval(r_vec, d)
  border <- ifelse(idx == 0, 0, cw[pmax(idx, 1)])

  return(data.frame(r = r_vec, border = 2 * border))
}


#' Use importance sampling to estimate the (stationary) intensity
#' @description
#' Using a thinned Thomas process simulated with parameters kappa_0, mu_0, and omega_0,
#' we use importance sampling to estimate
#'
#' @param w_is Importance sampling weights
#' @param patternSim List of point patterns to use for estimation
#' @param rho_baseline rho_baseline
#' @param normalized Set to TRUE if you want to normalize the importance sampling ratios
#'
#' @export
rho_importance_sampling <- function(w_is, rho_baseline, patternSim = NULL, normalized = FALSE){
  if(is.null(rho_baseline)){
    rho_baseline <- sapply(patternSim, estimate_rho_baseline)
  }

  if(normalized){
    w_is <- w_is/sum(w_is)
    return(sum(w_is*rho_baseline))
  } else {
    return(mean(w_is*rho_baseline))
  }
}

#' Estimate the intensity directly for a given point pattern
#' @description
#' Estimate the intensity directly for a given point pattern.
#' This estimate assumes stationarity.
#'
#' @param pattern The point pattern to estimate for
#'
#' @export
estimate_rho_baseline <- function(pattern){
  return(nrow(pattern$thinned)/
           ((pattern$xlim_thinned[2] - pattern$xlim_thinned[1])*(pattern$ylim_thinned[2] - pattern$ylim_thinned[1])))
}

#' Calculate Poisson density of a parent process
#' @description
#' Computes the log-density of a daughter point pattern \eqn{C} as a homogeneous
#' Poisson Point pattern restricted to the window \eqn{B}, with intensity \eqn{\kappa}.
#'
#' Specifically, the density is given by
#' \deqn{f(C\cup B|\kappa)=\exp\left(|B|(1-\kappa)\right)\kappa^{|C\cap B|}}
#'
#'
#' @param kappa    Positive scalar. Intensity of parent points.
#' @param parent   A \code{data.frame} with columns \code{x} and \code{y} giving parent locations,
#' typically restricted to an enlarged window \eqn{B \supset W}.
#' @param B_area   Area of the enlarged window \eqn{B \supset W}.
#' @param pattern  A list containing elements "parent" and "B_area" as described above
#'
#' @export
parent_log_density <- function(pattern = NULL, kappa, parent = NULL, B_area = NULL){
  return(pattern$B_area*(1-kappa) + nrow(pattern$parent)*log(kappa))
  # if(!is.null(pattern)){
  #   return(pattern$B_area*(1-kappa) + nrow(pattern$parent)*log(kappa))
  # } else {
  #   return(B_area*(1-kappa) + nrow(parent)*log(kappa))
  # }
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
  return(thomas_daughter_log_density_from_sums(mu = mu, kernel_sums = kernel_sums))
}

#' Kernel sums for the Thomas daughter log-density (the omega-dependent part)
#'
#' @description
#' Computes the parts of \code{\link{thomas_daughter_log_density}} that depend on
#' \eqn{\omega} (and the fixed parent/daughter positions) but not on \eqn{\mu}:
#' \eqn{\sum_j p_W(c_j)} and \eqn{\sum_i \log\left(\sum_j k(x_i-c_j)\right)}. These are
#' the expensive, \eqn{O(n_{daughter} \times n_{parent})} parts of the log-density. Split
#' out so they can be computed once and reused across evaluations that only change
#' \eqn{\mu}, via \code{\link{thomas_daughter_log_density_from_sums}}.
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
    K <- stats::dnorm(dx, mean = 0, sd = omega) * stats::dnorm(dy, mean = 0, sd = omega)
    sumLogK <- sum(log(rowSums(K)))
  }

  return(list(area_W = area_W, S_pW = S_pW, sumLogK = sumLogK, nDaughter = nrow(daughter)))
}

#' Combine cached kernel sums with mu into the Thomas daughter log-density
#'
#' @description
#' O(1) combine step completing \code{\link{thomas_daughter_kernel_sums}} into the
#' full log-density from \code{\link{thomas_daughter_log_density}}.
#'
#' @param mu Positive scalar. Mean number of daughters per parent.
#' @param kernel_sums Output of \code{\link{thomas_daughter_kernel_sums}}.
#'
#' @return A scalar giving \eqn{\log f(X \mid C)}.
#'
#' @export
thomas_daughter_log_density_from_sums <- function(mu, kernel_sums){
  kernel_sums$area_W - mu*kernel_sums$S_pW + kernel_sums$sumLogK + log(mu)*kernel_sums$nDaughter
}
