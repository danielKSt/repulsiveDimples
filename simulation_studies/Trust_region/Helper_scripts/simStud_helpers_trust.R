# Load libraries
suppressPackageStartupMessages({
  library(spatstat)
  library(ggplot2)
  library(parallel)
  library(mcprogress)
  library(repulsiveDimples)
})

generate_data <- function(sidelength, sim_pars, r_vec){
  pattern <- rThomas_matern_thinned(kappa = sim_pars$kappa,
                                    scale = sim_pars$omega,
                                    mu = sim_pars$mu,
                                    repulsionRange = sim_pars$rRange,
                                    xlims = c(0, sidelength),
                                    ylims = c(0, sidelength),
                                    saveparents = TRUE)
  rho_est <- estimate_rho_baseline(pattern = pattern)
  K_est <- estimate_K_lambda_baseline(pattern = pattern, r_vec = r_vec)
  return(list(rho = rho_est, K = K_est))
}
