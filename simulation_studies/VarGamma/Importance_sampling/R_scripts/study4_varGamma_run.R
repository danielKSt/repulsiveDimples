# Study 4 for the Variance gamma process: the design of study4_run.R in
# simulation_studies/Thomas/Importance_sampling/R_scripts, with the Thomas process
# replaced by a Variance gamma process with nu fixed at 1/2, where the kernel is the
# exponential kernel exp(-r/scale)/(2*pi*scale^2).
#
# The two par_goal groups keep study 4's kappa and mu, and their scale is set so that
# the clusters are the same size as there: at nu = 1/2 a daughter's mean squared
# distance from its parent is 6*scale^2, against 2*omega^2 for the Thomas process, so
# scale = omega/sqrt(3). The two studies then differ only in the shape of the kernel.
# The perturbation grids are study 4's, with scale on omega's coarser grid.
#
# This costs far more than study 4. The kernel's window probabilities need numerical
# integration, about 0.1 ms per parent, and the second group's patterns have about 370
# parents each. A full run is about 33 CPU-hours, most of it in the second group.
#
# The random number generator is set to L'Ecuyer-CMRG, since with the default one the
# forked mclapply workers reseed themselves and set.seed does not make the run
# reproducible. With it, the run reproduces given the same mc.cores.
#
# Run from the project root. The results are saved in
# simulation_studies/VarGamma/Importance_sampling/Results.
#
# For a fast end-to-end check, add --dry-run to the command line:
#   Rscript simulation_studies/VarGamma/Importance_sampling/R_scripts/study4_varGamma_run.R --dry-run

source("simulation_studies/VarGamma/Importance_sampling/R_scripts/importance_sampling_varGamma.R")

dry_run <- "--dry-run" %in% commandArgs(trailingOnly = TRUE)

repRange <- 0.4
nu <- 1/2

par_goal <- data.frame(kappa = log(c(0.32, 1.2)),
                       scale = log(c(0.8, 1.2)/sqrt(3)),
                       mu = log(c(5.2, 2.4)))

fine_grid  <- c(0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.8)
kappa_grid <- sort(c(-fine_grid, fine_grid))
mu_grid    <- sort(c(-fine_grid, fine_grid))

coarse_grid <- c(0.1, 0.2, 0.4, 0.8)
scale_grid  <- sort(c(-coarse_grid, coarse_grid))

par_delta <- rbind(
  data.frame(kappa = 0,          scale = 0,          mu = 0),          # shared delta=0 anchor
  data.frame(kappa = kappa_grid, scale = 0,          mu = 0),
  data.frame(kappa = 0,          scale = scale_grid, mu = 0),
  data.frame(kappa = 0,          scale = 0,          mu = mu_grid)
)

sidelength <- 4
xlims_sim  <- c(0, sidelength)
ylims_sim  <- c(0, sidelength)
xlims_true <- c(0, sidelength)
ylims_true <- c(0, sidelength)

r_vec <- seq(from = 0, to = sidelength * 3/4, by = 0.1)

if (dry_run) {
  nSims <- 200
  nTrue <- 200
} else {
  nSims <- 60000
  nTrue <- 200000
}

nSamples <- floor(10^((seq(from = log10(10), to = log10(nSims), by = 0.333334))))
if (nSamples[length(nSamples)] < nSims) { nSamples <- c(nSamples, nSims) }

cat(sprintf("study 4, Variance gamma (nu = %g): %d par_goal groups x %d delta cases = %d cases\n",
            nu, nrow(par_goal), nrow(par_delta), nrow(par_goal) * nrow(par_delta)))
cat(sprintf("  window %gx%g, nSims %d, nTrue %d, mc.cores %s%s\n",
            sidelength, sidelength, nSims, nTrue, getOption("mc.cores"),
            if (dry_run) "   [DRY RUN]" else ""))

RNGkind("L'Ecuyer-CMRG")
set.seed(16180)
simStudyResults <- combined_sim_study_varGamma(par_delta = par_delta, nu = nu,
                                               par_goal = par_goal, repRange = repRange,
                                               nSims = nSims, xlims_true = xlims_true, ylims_true = ylims_true,
                                               xlims_sim = xlims_sim, ylims_sim = ylims_sim, nSamples = nSamples,
                                               r_vec = r_vec, nTrue = nTrue,
                                               printProgress = TRUE)

resFolder <- "simulation_studies/VarGamma/Importance_sampling/Results/"
save(simStudyResults, file = paste0(resFolder, if (dry_run) "study4_varGamma_dryrun.RDa" else "study4_varGamma.RDa"))
