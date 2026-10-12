
# Study 3a: per-parameter perturbation curves.
# Same two par_goal groups as study2, but each of kappa/omega/mu is now
# perturbed one at a time across a grid of magnitudes (instead of just +-0.2),
# to get an actual dose-response curve of IS bias/degeneracy vs. shift size
# per parameter per regime. omega uses a coarser grid (it was very forgiving
# in study2). No parameters are shifted simultaneously here - that's left for
# a later interaction study.

source("importance_sampling_thomas.R")

repRange <- 0.4

par_goal <- data.frame(kappa = log(c(0.32, 1.2)),
                       omega = log(c(0.8, 1.2)),
                       mu = log(c(5.2, 2.4)))

kappa_grid <- seq(from = -2, to = 1.2, by = 0.4)
mu_grid    <- seq(from = -2, to = 1.2, by = 0.4)
omega_grid <- c(-2, -1, -0.5, 0.5, 1, 2)

kappa_grid <- kappa_grid[! kappa_grid == 0]
mu_grid <- mu_grid[! mu_grid == 0]

par_delta <- rbind(
  data.frame(kappa = 0,          omega = 0,          mu = 0),          # shared delta=0 anchor
  data.frame(kappa = kappa_grid, omega = 0,          mu = 0),
  data.frame(kappa = 0,          omega = omega_grid, mu = 0),
  data.frame(kappa = 0,          omega = 0,          mu = mu_grid)
)

nSims <- 60
xlims_sim <- c(0, 6)
ylims_sim <- c(0, 6)
nSamples <- floor(10^((seq(from = log10(10), to = log10(nSims), by = 0.333334))))
if (nSamples[length(nSamples)] < nSims) { nSamples <- c(nSamples, nSims) }
r_vec <- seq(from = 0, to = 3, by = 0.1)

nTrue <- 200
xlims_true <- c(0, 6)
ylims_true <- c(0, 6)

set.seed(31415)
simStudyResults <- combined_sim_study_thomas(par_delta = par_delta,
                                             par_goal = par_goal, repRange = repRange,
                                             nSims = nSims, xlims_true = xlims_true, ylims_true = ylims_true,
                                             xlims_sim = xlims_sim, ylims_sim = ylims_sim, nSamples = nSamples,
                                             r_vec = r_vec, nTrue = nTrue,
                                             printProgress = TRUE)

save(simStudyResults, file = "study3a.RDa")
