
source(file = "simulation_studies/Trust_region/simStud_helpers_trust.R")

# Parameters:
true_kappa <- 0.1234
true_scale <- 0.8
true_mu <- 2/true_kappa
true_rRange <- 0.1
par_thomas <- data.frame(kappa = true_kappa,
                         omega = true_scale,
                         mu = true_mu,
                         rRange = true_rRange)
rm(true_scale, true_kappa, true_mu, true_rRange)

# Simulate true data: ----
sidelength <- 40

nTrue <- 2000000
set.seed(seed = 1350)
sim_list_thomas <- mcprogress::pmclapply(X = rep(sidelength, nTrue), FUN = generate_data,
                                         sim_pars = par_thomas, r_vec = seq(from = 0, to = 3, by = 0.05),
                                         mc.cores = 6)

rho_baseline <- sapply(sim_list_thomas, function(res) res$rho)
rho_hat <- mean(rho_baseline)
# sd(rho_baseline)
# mean(rho_baseline)

K_lambda_baseline <- lapply(sim_list_thomas, function(res) res$K)
K_hat <- K_importance_sampling(w_is = rep(1, nTrue), K_lambda_baseline = K_lambda_baseline,
                               rho_baseline = rho_baseline)

plot(K_hat, type = "l")

# Want to also test with union approach:
nUnions <- 1000
sim_list_thomas <- mcprogress::pmclapply(X = rep(par_thomas$kappa, nUnions), FUN = rThomas_matern_thinned,
                                         scale = par_thomas$omega,
                                         mu = par_thomas$mu,
                                         repulsionRange = par_thomas$rRange,
                                         xlims = c(0, sidelength),
                                         ylims = c(0, sidelength),
                                         saveparents = TRUE,
                                         mc.cores = 6)
points_input <- lapply(X = sim_list_thomas, function(a) {
  if(nrow(a$thinned) == 0){
    return(0)
  } else {
    return(a$thinned)
  }
})
K_hat_unions <- K_est.unions(points_input = points_input, l = sidelength, spacing = 5,
                             r_vec = seq(from = 0, to = 3, by = 0.05), timescale = 1)

plot(x = K_hat$r, y = K_hat$border^(1/4), type = "l")
lines(x = K_hat_unions$r, y = K_hat_unions$border^(1/4), col = "red", lty = 2)

save(par_thomas, rho_hat, K_hat, K_hat_unions, nTrue, nUnions, sidelength, file = "study2.RDa")
