#' Use importance sampling to estimate the (stationary) intensity
#' @description
#' Using a thinned Thomas process simulated with parameters kappa_0, mu_0, and omega_0,
#' we use importance sampling to estimate
#'
#' @param w_is Importance sampling weights
#' @param patternSim List of point patterns to use for estimation. Only used when
#' `rho_baseline` is `NULL`, in which case the baseline intensities are estimated from
#' the patterns with \code{\link{estimate_rho_baseline}}.
#' @param rho_baseline Numeric vector of per-pattern baseline intensity estimates, one
#' per element of `patternSim`, as returned by \code{\link{estimate_rho_baseline}}. Pass
#' these in when they have already been computed (e.g. by
#' \code{\link{simulation_step}}) to avoid re-estimating them.
#' @param normalized Set to TRUE if you want to normalize the importance sampling ratios
#'
#' @return A scalar estimate of the stationary intensity under the target parameters.
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

#' Use importance sampling to estimate the K-function
#' @description
#' Use importance sampling to estimate the Ripleys' K-function for parameter values
#'
#' @param w_is Importance sampling weights
#' @param patternSim List of point patterns to use for estimation. Only used when
#' `K_lambda_baseline` or `rho_baseline` is `NULL`, in which case the missing baselines
#' are estimated from the patterns.
#' @param K_lambda_baseline List of per-pattern baseline estimates of the unnormalized
#' K-function \eqn{\lambda^2 K(r)}, one per element of `patternSim`, as returned by
#' \code{\link{estimate_K_lambda_baseline}}.
#' @param r_vec Vector of radii to estimate for. Only used when `K_lambda_baseline` is
#' `NULL` and the baselines have to be estimated from `patternSim`.
#' @param rho_baseline Numeric vector of per-pattern baseline intensity estimates, one
#' per element of `patternSim`, as returned by \code{\link{estimate_rho_baseline}}. Used
#' to divide the unnormalized K-function by the squared intensity.
#' @param normalized Set to TRUE if you want to normalize the importance sampling ratios
#'
#' @return A `data.frame` with columns `r` and `border`, giving the estimated K-function
#' under the target parameters.
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
#' @param patternSim List of point patterns to use for estimation. Only used when
#' `K_lambda_baseline` is `NULL`, in which case the baselines are estimated from the
#' patterns with \code{\link{estimate_K_lambda_baseline}}.
#' @param K_lambda_baseline List of per-pattern baseline estimates of the unnormalized
#' K-function \eqn{\lambda^2 K(r)}, one per element of `patternSim`, as returned by
#' \code{\link{estimate_K_lambda_baseline}}. Pass these in when they have already been
#' computed (e.g. by \code{\link{simulation_step}}) to avoid re-estimating them.
#' @param r_vec Vector of radii to estimate for. Only used when `K_lambda_baseline` is
#' `NULL` and the baselines have to be estimated from `patternSim`.
#' @param normalized Set to TRUE if you want to normalize the importance sampling ratios
#'
#' @return A `data.frame` with columns `r` and `border`, giving the estimated
#' unnormalized K-function under the target parameters.
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

#' Calculate Poisson density of a parent process
#' @description
#' Computes the log-density of a parent point pattern \eqn{C} as a homogeneous
#' Poisson point process restricted to the window \eqn{B}, with intensity \eqn{\kappa}.
#' The parents of every Neyman-Scott process are such a process, so this does not depend
#' on the cluster kernel.
#'
#' Specifically, the density is given by
#' \deqn{f(C\cap B|\kappa)=\exp\left(|B|(1-\kappa)\right)\kappa^{|C\cap B|}}
#'
#'
#' @param kappa    Positive scalar. Intensity of parent points.
#' @param parent   A \code{data.frame} with columns \code{x} and \code{y} giving parent locations,
#' typically restricted to an enlarged window \eqn{B \supset W}.
#' @param B_area   Area of the enlarged window \eqn{B \supset W}.
#' @param pattern  A list containing elements "parent" and "B_area" as described above. If
#' given, `parent` and `B_area` are taken from it.
#'
#' @export
parent_log_density <- function(pattern = NULL, kappa, parent = NULL, B_area = NULL){
  if(!is.null(pattern)){
    parent <- pattern$parent
    B_area <- pattern$B_area
  }
  return(B_area*(1-kappa) + nrow(parent)*log(kappa))
}

#' Combine cached kernel sums with mu into the daughter log-density
#'
#' @description
#' O(1) combine step completing a kernel's kernel sums into the log-density of the
#' daughter pattern \eqn{X} given the parents \eqn{C}. Given \eqn{C}, the daughters of a
#' Neyman-Scott process with a Poisson(\eqn{\mu}) number of daughters per parent and
#' dispersal kernel \eqn{k} form a Poisson process on \eqn{W} with intensity
#' \eqn{\mu \sum_{c_j \in C} k(x - c_j)}, so
#' \deqn{
#'   \log f(X \mid C) = |W| - \mu \sum_{c_j \in C} p_W(c_j)
#'   + \sum_{x_i \in X} \log \sum_{c_j \in C} k(x_i - c_j) + |X| \log \mu
#' }
#' with \eqn{p_W(c_j) = \int_W k(x - c_j)\,\mathrm{d}x}. Only the kernel sums depend on the
#' kernel, so this step is shared by every kernel. See
#' \code{\link{thomas_daughter_log_density}} for the Thomas process.
#'
#' @param mu Positive scalar. Mean number of daughters per parent.
#' @param kernel_sums A list with \code{area_W} (\eqn{|W|}), \code{S_pW}
#' (\eqn{\sum_j p_W(c_j)}), \code{sumLogK} (\eqn{\sum_i \log \sum_j k(x_i-c_j)}) and
#' \code{nDaughter} (\eqn{|X|}), as returned by \code{\link{thomas_daughter_kernel_sums}}.
#'
#' @return A scalar giving \eqn{\log f(X \mid C)}.
#'
#' @export
daughter_log_density_from_sums <- function(mu, kernel_sums){
  kernel_sums$area_W - mu*kernel_sums$S_pW + kernel_sums$sumLogK + log(mu)*kernel_sums$nDaughter
}
