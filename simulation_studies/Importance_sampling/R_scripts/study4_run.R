
# Study 4: study3a redone with a smaller window and smaller perturbations.
#
# study3a used a 6x6 window with log-deltas out to +-2, and the effective sample
# size collapsed: ESS was a median of ~0.01% of nSims in both parameter groups
# (i.e. of 60000 realizations the weighted ensemble was worth fewer than ten),
# so most of that grid says little beyond "the weights have degenerated".
#
# Two changes here:
#   * 4x4 window instead of 6x6, which roughly halves the number of latent
#     points each weight is a product over.
#   * log-deltas capped at +-0.8, with the grid concentrated near zero
#     (0.05, 0.1, 0.2, 0.4, 0.8) so the informative region is resolved rather
#     than spending most cases out where ESS is already gone.
#
# One parameter is perturbed at a time; simultaneous shifts are left to a later
# interaction study. omega keeps a coarser grid than kappa/mu, as before, and is
# now also capped at +-0.8 - study3a already covers wide omega (out to +-2), so
# there is no need to repeat it here.
#
# Set dry_run <- TRUE for a fast end-to-end check. The values below are the ones
# meant for the real run, so this script stays a faithful record of how
# study4.RDa was produced (study3a_run.R drifted to dry-run values and no longer
# reproduces study3a.RDa).

source("importance_sampling_thomas.R")

dry_run <- FALSE

repRange <- 0.4

par_goal <- data.frame(kappa = log(c(0.32, 1.2)),
                       omega = log(c(0.8, 1.2)),
                       mu = log(c(5.2, 2.4)))

# Doubling from 0.05 out to the 0.8 cap, mirrored about zero: dense where the
# weights should still be usable, sparse out at the edge.
fine_grid  <- c(0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.8)
kappa_grid <- sort(c(-fine_grid, fine_grid))
mu_grid    <- sort(c(-fine_grid, fine_grid))

coarse_grid <- c(0.1, 0.2, 0.4, 0.8)
omega_grid  <- sort(c(-coarse_grid, coarse_grid))

par_delta <- rbind(
  data.frame(kappa = 0,          omega = 0,          mu = 0),          # shared delta=0 anchor
  data.frame(kappa = kappa_grid, omega = 0,          mu = 0),
  data.frame(kappa = 0,          omega = omega_grid, mu = 0),
  data.frame(kappa = 0,          omega = 0,          mu = mu_grid)
)

sidelength <- 4
xlims_sim  <- c(0, sidelength)
ylims_sim  <- c(0, sidelength)
xlims_true <- c(0, sidelength)
ylims_true <- c(0, sidelength)

# Keep rMax at 3/4 of the side, matching the rule used in the window study.
r_vec <- seq(from = 0, to = sidelength * 3/4, by = 0.1)

if (dry_run) {
  nSims <- 200
  nTrue <- 200
} else {
  nSims <- 60000    # as in study3a, so ESS is comparable at equal sample size
  nTrue <- 200000   # smaller window makes this cheaper than study3a's 200000 at 6x6
}

nSamples <- floor(10^((seq(from = log10(10), to = log10(nSims), by = 0.333334))))
if (nSamples[length(nSamples)] < nSims) { nSamples <- c(nSamples, nSims) }

cat(sprintf("study 4: %d par_goal groups x %d delta cases = %d cases\n",
            nrow(par_goal), nrow(par_delta), nrow(par_goal) * nrow(par_delta)))
cat(sprintf("  window %gx%g, nSims %d, nTrue %d, mc.cores %s%s\n",
            sidelength, sidelength, nSims, nTrue, getOption("mc.cores"),
            if (dry_run) "   [DRY RUN]" else ""))

set.seed(27182)
simStudyResults <- combined_sim_study_thomas(par_delta = par_delta,
                                             par_goal = par_goal, repRange = repRange,
                                             nSims = nSims, xlims_true = xlims_true, ylims_true = ylims_true,
                                             xlims_sim = xlims_sim, ylims_sim = ylims_sim, nSamples = nSamples,
                                             r_vec = r_vec, nTrue = nTrue,
                                             printProgress = TRUE)

save(simStudyResults, file = if (dry_run) "study4_dryrun.RDa" else "study4.RDa")
