#' Objective function for minimum contrast
#' @description
#' Calculate contrast/objective function given params and estimates for the intensity and the K-function from data
#'
#'
#' @param params Parameter vector(to be optimized)
#' @param rho_hat Intensity estimated from data
#' @param K_hat K-function estimated from data
#' @param params_fixed Values for the fixed parameters
#' @param free_params_indices Which indices for the parameter vector are free?
#' @param win Simulation window
#' @param nSims How many simulations to estimate K_hat
#' @param wpq Vector of weights and powers (w_rho, q_rho, p_rho, w_K, q_K, p_K)
#' @param ns_kernel String indicating which kernel to use for the Neymann-Scott process, "var_gamma" or "thomas"
#' @param union_est Set to TRUE in order to use union approach for estimation of K-function
#' @param parallel Set to TRUE to use parallel computing when simulating point patterns for estimation
#'
#' @export
objective_weighted_contrast <- function(params, rho_hat, K_hat, params_fixed, free_params_indices = NULL,
                                                 wpq = c(1, 1, 2, 1, 1/4, 2),
                                                 ns_kernel = "var_gamma",
                                                 nSims, win, union_est = FALSE,
                                                 parallel = TRUE){
  return(0)
  # if(!is.null(free_params_indices)){
  #   params_temp <- c(1:(length(params)+length(params_fixed)))
  #   params_temp[free_params_indices] <- params
  #   params <- params_temp
  #   params[-free_params_indices] <- params_fixed
  #   rm(params_temp)
  # }
  #
  # # if(is.numeric(repulsionRange_hat)){
  # #   params <- c(params, repulsionRange_hat)
  # # }
  # r_vec <- K_hat$r
  # if(ns_kernel == "var_gamma"){
  #   simulated_summaries <- K_theo.var_gamma(r_vec = r_vec, params = params, nSims = nSims,
  #                                           win = win, union_est = union_est, parallel = parallel)
  # } else if(ns_kernel == "thomas"){
  #   simulated_summaries <- K_theo.var_gamma(r_vec = r_vec, params = params, nSims = nSims,
  #                                           win = win, union_est = union_est, parallel = parallel)
  # } else {
  #   simulated_summaries <- K_theo.var_gamma(r_vec = r_vec, params = params, nSims = nSims,
  #                                           win = win, union_est = union_est, parallel = parallel)
  # }
  #
  #
  # rho_theo <- simulated_summaries$rho_hat
  # if(union_est){
  #   K_theo_est <- simulated_summaries$K_hat$border
  # } else {
  #   K_theo_est <- simulated_summaries$K_hat$poolborder
  # }
  #
  # dr <- r_vec[2] - r_vec[1]
  # contrast_K <- sum(abs(K_hat$border[which(!is.na(K_theo_est))]^wpq[6] - K_theo_est[which(!is.na(K_theo_est))]^wpq[6])^wpq[5])*dr
  # contrast_rho <- abs(rho_theo^wpq[3] - rho_hat^wpq[3])^wpq[2]
  #
  # return(wpq[1]*contrast_rho + wpq[4]*contrast_K)
}

#' Estimate K-function for variance gamma point pattern
#' @description
#' To be written
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
K_theo.var_gamma <- function(r_vec, params, nSims, win, union_est = FALSE, parallel = TRUE, nCores = 4){
  return(0)
  # varphi <- exp(params[1])
  # rho <- exp(params[2])
  # sigmasq <- exp(params[3])
  # kappa <- rho^2/(2*pi*varphi*sigmasq)
  # mu <- 2*pi*varphi*sigmasq/rho
  # rRange <- exp(params[4])
  # if(parallel){
  #   sims <- parallel::mclapply(X = rep(kappa, nSims),
  #                              FUN = function(kappa, scale, mu, nu, rRange, win){
  #                                a <-rVarGamma_matern_thinned(kappa = kappa, scale = varphi, mu = mu, nu = -1/4,
  #                                                             repulsionRange = rRange, win = win)
  #                                if(a$n > 0){
  #                                  return(data.frame(x = a$x/win$xrange[2], y = a$y/win$yrange[2]))
  #                                } else {
  #                                  return(0)
  #                                }
  #                                },
  #                              scale = varphi,
  #                              mu = mu,
  #                              nu = -1/4,
  #                              rRange = rRange,
  #                              win = win,
  #                              mc.cores = nCores)
  # } else {
  #   my_expr <- function(){
  #       a <-rVarGamma_matern_thinned(kappa = kappa, scale = varphi, mu = mu, nu = -1/4,
  #                                    repulsionRange = rRange, win = win)
  #       if(a$n > 0){
  #           return(data.frame(x = a$x/win$xrange[2], y = a$y/win$yrange[2]))
  #       } else {
  #         return(0)
  #       }
  #     }
  #   sims <- replicate(n = nSims,
  #                     expr = my_expr(),
  #                     simplify = FALSE)
  # }
  # nTot <- sum(sapply(sims, FUN = function(x){
  #   if(is.numeric(x)){return(0)}
  #   else{return(nrow(x))}}))
  # if(union_est){
  #   return(list(K_hat = K_est.unions(points_input = sims,
  #                                    l = win$xrange[2],
  #                                    spacing = 5,
  #                                    r_vec = r_vec,
  #                                    timescale = 1),
  #               rho_hat = nTot/(spatstat.geom::area(win)*nSims)))
  # } else {
  #   Keach <- lapply(sims, spatstat.explore::Kest, r = r_vec)
  #   return(list(K_hat = spatstat.explore::pool(spatstat.geom::as.anylist(Keach)),
  #               rho_hat = nTot/(spatstat.geom::area(win)*nSims)))
  # }
}

#' Estimate K-function for thinned Thomas process
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
K_theo.thomas <- function(r_vec, params, nSims, win, union_est = FALSE, parallel = TRUE, nCores = 4){
  return(0)
  # mom_rho <- exp(params[1])
  # omega <- exp(params[2])
  # mu <- exp(params[3])
  # rRange <- exp(params[4])
  # if(parallel){
  #   sims <- parallel::mclapply(X = rep(kappa, nSims),
  #                              FUN = function(mom_rho, omega, mu, rRange, win){
  #                                a <-rThomas_matern_thinned(kappa = mom_rho, scale = omega,
  #                                                           mu = mu, repulsionRange = rRange, win = win)
  #                                if(a$n > 0){
  #                                  return(data.frame(x = a$x/win$xrange[2], y = a$y/win$yrange[2]))
  #                                } else {
  #                                  return(0)
  #                                }
  #                              },
  #                              scale = omega,
  #                              mu = mu,
  #                              repulsionRange = rRange,
  #                              win = win,
  #                              mc.cores = nCores)
  # } else {
  #   my_expr <- function(){
  #     a <-rThomas_matern_thinned(kappa = mom_rho, scale = omega, mu = mu,
  #                                repulsionRange = rRange, win = win)
  #     if(a$n > 0){
  #       return(data.frame(x = a$x/win$xrange[2], y = a$y/win$yrange[2]))
  #     } else {
  #       return(0)
  #     }
  #   }
  #   sims <- replicate(n = nSims,
  #                     expr = my_expr(),
  #                     simplify = FALSE)
  # }
  # nTot <- sum(sapply(sims, FUN = function(x){
  #   if(is.numeric(x)){return(0)}
  #   else{return(nrow(x))}}))
  # if(union_est){
  #   return(list(K_hat = K_est.unions(points_input = sims,
  #                                    l = win$xrange[2],
  #                                    spacing = 5,
  #                                    r_vec = r_vec,
  #                                    timescale = 1),
  #               rho_hat = nTot/(spatstat.geom::area(win)*nSims)))
  # } else {
  #   Keach <- lapply(sims, spatstat.explore::Kest, r = r_vec)
  #   return(list(K_hat = spatstat.explore::pool(spatstat.geom::as.anylist(Keach)),
  #               rho_hat = nTot/(spatstat.geom::area(win)*nSims)))
  # }
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
