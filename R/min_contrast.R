# Variance gamma contrast fitting ----

#' Minimal contrast for unconstrained estimation with thinned varGamma
#' @description
#' Use minimal contrast estimation together with Monte Carlo calculations for K_theo for parameter estimation for Matern thinned VarGamma
#'
#'
#' @param params_init Initial guess for parameter vector
#' @param rho_hat Intensity estimated from data
#' @param K_hat K-function estimated from data
#' @param repulsionRange_hat Set to FALSE if the repulsion range parameter is in the params vector
#' @param win Simulation window
#' @param nSims How many simulations to estimate K_hat
#' @param union_est Set to TRUE in order to use union approach for estimation of K-function
#' @param w_rho Weight of intensity
#' @param q_rho Power to use inside parantheses in contrast for intensity
#' @param p_rho Power to use outside parantheses in contrast for intensity
#' @param w_K Weight of K-function
#' @param q_K Power to use inside parantheses in contrast for K-function
#' @param p_K Power to use outside parantheses in contrast for K-function
#'
#' @export
contrast_estimation_vg.unconstrained <- function(params_init, rho_hat, K_hat, repulsionRange_hat,
                                                 nSims, win, union_est = FALSE,
                                                 w_rho = 1, q_rho = 1, p_rho = 2,
                                                 w_K = 2, q_K = 1/4, p_K = 2){
  res <- stats::optim(par = params_init,
               fn = objective_weighted_contrast_vargamma,
               rho_hat = rho_hat,
               K_hat = K_hat,
               repulsionRange_hat = repulsionRange_hat,
               w_rho = w_rho,
               q_rho = q_rho,
               p_rho = p_rho,
               w_K = w_K,
               q_K = q_K,
               p_K = p_K,
               union_est = union_est,
               win = win,
               nSims = nSims)
  return(res)
}

#' Minimal contrast for constrained estimation with thinned varGamma
#' @description
#' Use minimal contrast estimation together with Monte Carlo calculations for K_theo for parameter estimation for Matern thinned VarGamma, with linear constraints on the parameters.
#' The parameter vector should be of the form (clusterRange, unthinnedIntensity, sigmaSquared, repulsionRange)
#'
#'
#' @param params_init Initial guess for parameter vector
#' @param rho_hat Intensity estimated from data
#' @param K_hat K-function estimated from data
#' @param win Simulation window
#' @param nSims How many simulations to estimate K_hat
#' @param union_est Set to TRUE in order to use union approach for estimation of K-function
#' @param ui constraint matrix(see constrOptim documentation)
#' @param ci constraint vector(see constrOptim documentation)
#' @param w_rho Weight of intensity
#' @param q_rho Power to use inside parantheses in contrast for intensity
#' @param p_rho Power to use outside parantheses in contrast for intensity
#' @param w_K Weight of K-function
#' @param q_K Power to use inside parantheses in contrast for K-function
#' @param p_K Power to use outside parantheses in contrast for K-function
#'
#' @export
contrast_estimation_vg.constrained <- function(params_init, rho_hat, K_hat, nSims, win, union_est, ui, ci,
                                                 w_rho = 1, q_rho = 1, p_rho = 2, w_K = 2, q_K = 1/4, p_K = 2){
  res <- stats::constrOptim(theta = params_init,
                      f = objective_weighted_contrast_vargamma,
                      ui = ui,
                      ci = ci,
                      rho_hat = rho_hat,
                      K_hat = K_hat,
                      repulsionRange_hat = FALSE,
                      w_rho = w_rho,
                      q_rho = q_rho,
                      p_rho = p_rho,
                      w_K = w_K,
                      q_K = q_K,
                      p_K = p_K,
                      union_est = union_est,
                      win = win,
                      nSims = nSims)
  return(res)
}


#' Objective function for minimum contrast
#' @description
#' Calculate contrast/objective function given params and estimates for the intensity and the K-function from data
#'
#'
#' @param params Parameter vector(to be optimized)
#' @param rho_hat Intensity estimated from data
#' @param K_hat K-function estimated from data
#' @param repulsionRange_hat Set to FALSE if the repulsion range parameter is in the params vector
#' @param win Simulation window
#' @param nSims How many simulations to estimate K_hat
#' @param w_rho Weight of intensity
#' @param q_rho Power to use inside parantheses in contrast for intensity
#' @param p_rho Power to use outside parantheses in contrast for intensity
#' @param w_K Weight of K-function
#' @param q_K Power to use inside parantheses in contrast for K-function
#' @param p_K Power to use outside parantheses in contrast for K-function
#' @param union_est Set to TRUE in order to use union approach for estimation of K-function
#'
#' @export
objective_weighted_contrast_vargamma <- function(params, rho_hat, K_hat, repulsionRange_hat,
                                                 w_rho = 1, q_rho = 1, p_rho = 2,
                                                 w_K = 1, q_K = 1/4, p_K = 2,
                                                 nSims, win, union_est = FALSE){
  if(is.numeric(repulsionRange_hat)){
    params <- c(params, repulsionRange_hat)
  }
  r_vec <- K_hat$r
  simulated_summaries <- K_theo.exp(r_vec = r_vec, params = params, nSims = nSims,
                                    win = win, union_est = union_est)

  rho_theo <- simulated_summaries$rho_hat
  if(union_est){
    K_theo_est <- simulated_summaries$K_hat$border
  } else {
    K_theo_est <- simulated_summaries$K_hat$poolborder
  }

  dr <- r_vec[2] - r_vec[1]
  contrast_K <- sum(abs(K_hat$border[which(!is.na(K_theo_est))]^q_K - K_theo_est[which(!is.na(K_theo_est))]^q_K)^p_K)*dr
  contrast_rho <- abs(rho_theo^q_rho - rho_hat^q_rho)^p_rho

  return(w_rho*contrast_rho + w_K*contrast_K)
}

#' Estimate K-function with exponential kernel
#' @description
#' To be written
#'
#'
#' @param params Parameter vector(to be optimized) should contain (range, intensity, variance, repulsionRange)
#' @param r_vec Vector with radii to use for contrast calculation
#' @param nSims Number of simulations for estimation
#' @param win Simulation window
#' @param parallel Set to TRUE to run simulations in parallel
#' @param nCores Number of cores for parallel
#' @param union_est Set to TRUE in order to use union approach for estimation of K-function
#'
#' @export
K_theo.exp <- function(r_vec, params, nSims, win, union_est = FALSE, parallel = TRUE, nCores = 4){
  varphi <- params[1]
  rho <- params[2]
  sigmasq <- params[3]
  kappa <- rho^2/(2*pi*varphi*sigmasq)
  mu <- 2*pi*varphi*sigmasq/rho
  if(parallel){
    sims <- parallel::mclapply(X = rep(kappa, nSims),
                               FUN = rVarGamma_matern_thinned,
                               scale = varphi,
                               mu = mu,
                               nu = -1/4,
                               repulsionRange = params[4],
                               win = win,
                               mc.cores = nCores)
  } else {
    sims <- replicate(n = nSims, expr = rVarGamma_matern_thinned(kappa = kappa,
                                                                 scale = varphi,
                                                                 mu = mu,
                                                                 nu = -1/4,
                                                                 repulsionRange = params[4],
                                                                 win = win),
                      simplify = FALSE)
  }
  nTot <- sum(sapply(sims, FUN = function(x){return(x$n)}))
  if(union_est){
    return(list(K_hat = K_est.unions(points_input = sims,
                                     l = win$xrange[2],
                                     spacing = 5,
                                     r_vec = r_vec,
                                     timescale = 1),
                rho_hat = nTot/(spatstat.geom::area(win)*nSims)))
  } else {
    Keach <- lapply(sims, spatstat.explore::Kest, r = r_vec)
    return(list(K_hat = spatstat.explore::pool(spatstat.geom::as.anylist(Keach)),
                rho_hat = nTot/(spatstat.geom::area(win)*nSims)))
  }
}


# Regression model contrast fitting ----
#
# objective_function <- function(theta, lambda_dimp, pcf_dimp, r_low, r_up){
#   r_ind <- which((pcf_dimp$r > r_low)*(pcf_dimp$r < r_up) == 1)
#   radiar <- pcf_dimp$r[r_ind]
#   dr = radiar[2] - radiar[1]
#   numerator <- theta[2]^2*theta[4]^2*exp(-radiar/theta[5])
#   denominator <- theta[1]^2+2*theta[1]*theta[3]+theta[3]^2
#   g_r <- 1 + numerator/denominator
#   res <- (theta[1]+theta[3] - lambda_dimp)^2 + sum((pcf_dimp$trans[r_ind]-g_r)^2)*dr
#   return(res)
# }
#
# res <- optim(par = c(0.5, 2.5, 1, 2, 1), fn = objective_function, lambda_dimp = pcf.dimples$lHat, pcf_dimp = pcf.dimples$g_r, r_low = 0.5, r_up = 2.5)
# res$par
#
# rm(simulatedScars, simulatedVortices, res, pcf.dimples,pcf.scars, settings, lengthscale, timescale, objective_function)
# load_case <- function(settings){
#   load(file = paste("data/", settings, "/simulatedPoints.RDa", sep = ""))
#   load(file = paste("data/", settings, "/scales.RDa", sep = ""))
#   return(list(scars = simulatedScars, dimples = simulatedVortices, lengthscale = lengthscale, timescale = timescale))
# }
#
# get_pcf_case <- function(settings, params){
#   points_in <- load_case(settings = settings$settings)
#   pcf_scars <- get_pcf(points_input = points_in$scars,
#                        l = settings$size/points_in$lengthscale,
#                        spacing = params$spacingScars,
#                        bw = params$bwScars,
#                        rMax = params$rMaxScars,
#                        dr = params$drScars,
#                        timescale = points_in$timescale)
#   pcf_dimples <- get_pcf(points_input = points_in$dimples,
#                          l = settings$size/points_in$lengthscale,
#                          spacing = params$spacingDimples,
#                          bw = params$bwDimples,
#                          rMax = params$rMaxDimples,
#                          dr = params$drDimples,
#                          timescale = points_in$timescale)
#   return(list(res.dimples = pcf_dimples, res.scars = pcf_scars))
# }
#
# settings_list <- list(list(settings = paste("RE", 2500, "_WE", "inf", sep = ""), size = 256), list(settings = paste("RE", 1000, "_WE", "inf", sep = ""), size = 128),
#                       list(settings = paste("RE", 2500, "_WE", 10, sep = ""), size = 256), list(settings = paste("RE", 1000, "_WE", 10, sep = ""), size = 128))
# params <- list(spacingScars = 5,
#                bwScars = 0.08,
#                rMaxScars = 3,
#                drScars = 0.05,
#                spacingDimples = 5,
#                bwDimples = 0.04,
#                rMaxDimples = 3,
#                drDimples = 0.01)
# res_pcf <- lapply(settings_list, get_pcf_case,
#                   params = params)
#
# scar.pcf <- vector(mode = "list", length = length(settings_list))
# dimples.pcf <- vector(mode = "list", length = length(settings_list))
# scar.lHat <- vector(mode = "list", length = length(settings_list))
# dimples.lHat <- vector(mode = "list", length = length(settings_list))
#
# for(i in 1:length(settings_list)){
#   scar.pcf[[i]] <- res_pcf[[i]]$res.scars$g_r
#   scar.lHat[[i]] <- res_pcf[[i]]$res.scars$lHat
#   dimples.pcf[[i]] <- res_pcf[[i]]$res.dimples$g_r
#   dimples.lHat[[i]] <- res_pcf[[i]]$res.dimples$lHat
# }
#
# rm(params, res_pcf, i)
#
# objective_function_multicase <- function(theta, lHat, gHat, r_low, r_up){
#   res <- 0
#   for (case in 1:length(lHat)) {
#     r_ind <- which((gHat[[case]]$r > r_low[case])*(gHat[[case]]$r < r_up[case]) == 1)
#     radiar <- gHat[[case]]$r[r_ind]
#     numerator <- exp(-radiar/theta[5+3*(case-1)])*(theta[2]*theta[4+3*(case-1)])^2
#     denominator <- theta[1]^2+2*theta[1]*theta[3+3*(case-1)]+theta[2+3*(case-1)]^2
#     g_r <- 1 + numerator/denominator
#     res <- res + (theta[1]+theta[3+3*(case-1)] - lHat[[case]])^2 + sum((gHat[[case]]$trans[r_ind]-g_r)^2)
#   }
#   return(res)
# }
#
# res <- optim(par = c(1, 1, rep(1, 3*length(dimples.lHat))), fn = objective_function_multicase,
#              lHat = dimples.lHat,
#              gHat = dimples.pcf,
#              r_low = list(0.5, 0.5, 0.5, 0.5),
#              r_up = list(3, 2, 3, 2))
#
# coeffs <- res$par[c(1,2)]
# df <- data.frame(Re = c(2500, 1000, 2500, 1000), We = c(999999, 999999, 10, 10),
#                  phiHat = res$par[c(3, 6, 9, 12)],
#                  sigmaHat = res$par[c(4, 7, 10, 13)],
#                  rhoHat = res$par[c(5, 8, 11, 14)])
