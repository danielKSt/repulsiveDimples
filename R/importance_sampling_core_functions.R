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
