
#' Contrast function
#' @description
#' Calculate contrast given params and an estimated K-function from data
#'
#'
#' @param params Parameter vector(to be optimized)
#' @param K_est K-function estimated from data
#' @param K_theo Function for how to calculate K_theo
#' @param r_vec Vector with radii to use for contrast calculation
#' @param dr Step size for r_vec
#' @param win Simulation window
#' @param q Power to use in contrast(inner)
#' @param p Power to use in contrast(outer)
#'
#' @export
contrast_theo_est <- function(params, K_est, K_theo, r_vec, dr, win, q = 1/4, p = 2){
  valueK_theo <- K_theo(r_vec, params, win = win)
  res <- sum(abs(valueK_theo^q-K_est^q)^p)*dr
  return(res)
}

#' Estimate K-function with exponential kernel
#' @description
#' To be written
#'
#'
#' @param params Parameter vector(to be optimized)
#' @param r_vec Vector with radii to use for contrast calculation
#' @param win Simulation window
#'
#' @export
K_theo.exp <- function(r_vec, params, win){
  varphi <- params$range
  kappa <- params$rho^2/(2*pi*varphi*params$sigmasq)
  mu <- 2*pi*varphi*params$sigmasq/params$rho
  sims <- replicate(n = params$nSims, expr = rVarGamma_matern_thinned(kappa = kappa,
                                                                      scale = range,
                                                                      mu = mu,
                                                                      nu = -1/4,
                                                                      repulsionRange = params$repulisionRange,
                                                                      win = win),
                    simplify = FALSE)
  Keach <- lapply(sims, spatstat.explore::Kest, r = r_vec)
  return(spatstat.explore::pool(spatstat.geom::as.anylist(Keach))$poolborder)
}
