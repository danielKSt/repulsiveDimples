# Perturbing nu for the Variance gamma process: the design of study4_varGamma_run.R,
# with the same two groups of (kappa, scale, mu) and the same window, but with the
# shape nu perturbed instead, around nu = 1/2.
#
# nu is perturbed additively, on its own scale rather than the log scale, since it can
# be negative. The kernel needs nu > -1, but spatstat's rVarGamma refuses nu <= -1/2,
# because the default padding of its parent window comes from clusterradius(), which
# needs it. So the simulation nu_0 = 1/2 + delta needs delta > -1. At fixed scale, nu
# changes the size of the clusters as well as the shape of the kernel: a daughter's mean
# squared distance from its parent is 4*scale^2*(nu + 1).
#
# The grid is asymmetric. A pilot with 3000 patterns per case put the ESS at 25% and
# 51% (groups 1 and 2) at delta = -0.4 but 0.25% and 2.2% at -0.8, while at +1.2 it was
# still 6.6% and 41%. So the negative side stops at -0.9, and the positive side runs
# out to +1.6.
#
# Every case needs the kernel sums at both nu_0 and nu, and the patterns at larger nu_0
# have more parents, so this costs a little more than study 4: about 5 hours on 10 cores.
#
# The random number generator is set to L'Ecuyer-CMRG, since with the default one the
# forked mclapply workers reseed themselves and set.seed does not make the run
# reproducible. With it, the run reproduces given the same mc.cores.
#
# Run from the project root. The results are saved in
# simulation_studies/VarGamma/Importance_sampling/Results.
#
# For a fast end-to-end check, add --dry-run to the command line:
#   Rscript simulation_studies/VarGamma/Importance_sampling/R_scripts/studyNu_varGamma_run.R --dry-run

source("simulation_studies/VarGamma/Importance_sampling/R_scripts/importance_sampling_varGamma.R")

dry_run <- "--dry-run" %in% commandArgs(trailingOnly = TRUE)

repRange <- 0.4

par_goal <- data.frame(kappa = log(c(0.32, 1.2)),
                       scale = log(c(0.8, 1.2)/sqrt(3)),
                       mu = log(c(5.2, 2.4)),
                       nu = 1/2)

nu_grid <- c(-0.9, -0.8, -0.6, -0.4, -0.3, -0.2, -0.15, -0.1, -0.05,
             0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.2, 1.6)

par_delta <- data.frame(kappa = 0, scale = 0, mu = 0, nu = c(0, nu_grid))   # first row: shared delta=0 anchor

if (any(outer(par_goal$nu, par_delta$nu, "+") <= -1/2)) {
  stop("Every simulation nu must be greater than -1/2, which rVarGamma requires.")
}

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

cat(sprintf("nu study, Variance gamma (nu = %g): %d par_goal groups x %d delta cases = %d cases\n",
            par_goal$nu[1], nrow(par_goal), nrow(par_delta), nrow(par_goal) * nrow(par_delta)))
cat(sprintf("  window %gx%g, nSims %d, nTrue %d, mc.cores %s%s\n",
            sidelength, sidelength, nSims, nTrue, getOption("mc.cores"),
            if (dry_run) "   [DRY RUN]" else ""))

RNGkind("L'Ecuyer-CMRG")
set.seed(57721)
simStudyResults <- combined_sim_study_varGamma(par_delta = par_delta,
                                               par_goal = par_goal, repRange = repRange,
                                               nSims = nSims, xlims_true = xlims_true, ylims_true = ylims_true,
                                               xlims_sim = xlims_sim, ylims_sim = ylims_sim, nSamples = nSamples,
                                               r_vec = r_vec, nTrue = nTrue,
                                               printProgress = TRUE)

resFolder <- "simulation_studies/VarGamma/Importance_sampling/Results/"
save(simStudyResults, file = paste0(resFolder, if (dry_run) "studyNu_varGamma_dryrun.RDa" else "studyNu_varGamma.RDa"))
