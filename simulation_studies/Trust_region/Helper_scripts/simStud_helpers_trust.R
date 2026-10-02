# Load libraries
suppressPackageStartupMessages({
  library(spatstat)
  library(ggplot2)
  library(parallel)
  library(mcprogress)
  library(repulsiveDimples)
})


generate_data <- function(sidelength, sim_pars, r_vec, thinningType = 2){
  pattern <- rThomas_matern_thinned(kappa = sim_pars$kappa,
                                    scale = sim_pars$omega,
                                    mu = sim_pars$mu,
                                    repulsionRange = sim_pars$rRange,
                                    xlims = c(0, sidelength),
                                    ylims = c(0, sidelength),
                                    thinningType = thinningType,
                                    saveparents = TRUE)
  rho_est <- estimate_rho_baseline(pattern = pattern)
  K_est <- estimate_K_lambda_baseline(pattern = pattern, r_vec = r_vec)
  return(list(rho = rho_est, K = K_est))
}

simulate_for_study_settings <- function(parThomas, studyNr, nTrue, nUnions, sidelength, thinningType = 2,
                                        mc.cores = 6, seed = 1350){
  set.seed(seed = seed)
  sim_list_thomas <- mcprogress::pmclapply(X = rep(sidelength, nTrue), FUN = generate_data,
                                           sim_pars = parThomas, r_vec = seq(from = 0, to = 3, by = 0.05),
                                           thinningType = thinningType, mc.cores = mc.cores)

  rho_baseline <- sapply(sim_list_thomas, function(res) res$rho)
  rho_hat <- mean(rho_baseline)

  K_lambda_baseline <- lapply(sim_list_thomas, function(res) res$K)
  K_hat <- K_importance_sampling(w_is = rep(1, nTrue), K_lambda_baseline = K_lambda_baseline,
                                 rho_baseline = rho_baseline)

  # Want to also test with union approach:
  sim_list_thomas <- mcprogress::pmclapply(X = rep(parThomas$kappa, nUnions), FUN = rThomas_matern_thinned,
                                           scale = parThomas$omega,
                                           mu = parThomas$mu,
                                           repulsionRange = parThomas$rRange,
                                           xlims = c(0, sidelength),
                                           ylims = c(0, sidelength),
                                           saveparents = TRUE,
                                           thinningType = thinningType,
                                           mc.cores = mc.cores)
  points_input <- lapply(X = sim_list_thomas, function(a) {
    if(nrow(a$thinned) == 0){
      return(0)
    } else {
      return(a$thinned)
    }
  })
  K_hat_unions <- K_est.unions(points_input = points_input, l = sidelength, spacing = 5,
                               r_vec = seq(from = 0, to = 3, by = 0.05), timescale = 1)


  # Saved as par_thomas, the name run_study and the existing data files use.
  par_thomas <- parThomas
  save(par_thomas, rho_hat, K_hat, K_hat_unions, nTrue, nUnions, sidelength,
       file = paste0("simulation_studies/Trust_region/Thinning_type", thinningType, "/Data/study", studyNr, ".RDa"))
}
