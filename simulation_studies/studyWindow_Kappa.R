
source("importance_sampling_thomas.R")

# Varying kappa: ----
repRange <- 0.4
sidelength <- 6
xlims_sim <- c(0, sidelength)
ylims_sim <- c(0, sidelength)
xlims_true <- c(0, sidelength)
ylims_true <- c(0, sidelength)
r_vec <- seq(from = 0, to = 3, by = 0.1)

nSims <- 60000
nSamples <- floor(10^((seq(from = log10(10), to = log10(nSims), by = 0.333334))))
if (nSamples[length(nSamples)] < nSims) { nSamples <- c(nSamples, nSims) }
nTrue <- 20000

par_goal <- data.frame(kappa = log(seq(from = 0.2, to = 2.0, by = 0.15)),
                       omega = log(c(0.8)),
                       mu = log(c(4.2)))

par_delta <- rbind(
  data.frame(kappa = c(-0.5, 0.5), omega = 0, mu = 0),
  data.frame(kappa = 0, omega = 0, mu = c(-0.5, 0.5))
)

set.seed(31415)
simStudyResults <- combined_sim_study_thomas(par_delta = par_delta,
                                             par_goal = par_goal, repRange = repRange,
                                             nSims = nSims, xlims_true = xlims_true, ylims_true = ylims_true,
                                             xlims_sim = xlims_sim, ylims_sim = ylims_sim, nSamples = nSamples,
                                             r_vec = r_vec, nTrue = nTrue,
                                             printProgress = TRUE)

save(simStudyResults, file = "studyKappa.RDa")

# Varying window sizes: ----

repRange <- 0.4

par_goal <- data.frame(kappa = log(c(0.32)),
                       omega = log(c(0.8)),
                       mu = log(c(4.2)))

par_delta <- rbind(
  data.frame(kappa = c(-0.5, 0.5), omega = 0, mu = 0),
  data.frame(kappa = 0, omega = 0, mu = c(-0.5, 0.5))
)

nSims <- 60000
nSamples <- floor(10^((seq(from = log10(10), to = log10(nSims), by = 0.333334))))
if (nSamples[length(nSamples)] < nSims) { nSamples <- c(nSamples, nSims) }
nTrue <- 20000

analysis_given_sidelength <- function(sidelength, rmax){
  r_vec <- seq(from = 0, to = min(rmax, sidelength*3/4), by = 0.1)

  xlims_sim <- c(0, sidelength)
  ylims_sim <- c(0, sidelength)
  xlims_true <- c(0, sidelength)
  ylims_true <- c(0, sidelength)

  simStudyResults <- combined_sim_study_thomas(par_delta = par_delta,
                                               par_goal = par_goal, repRange = repRange,
                                               nSims = nSims, xlims_true = xlims_true, ylims_true = ylims_true,
                                               xlims_sim = xlims_sim, ylims_sim = ylims_sim, nSamples = nSamples,
                                               r_vec = r_vec, nTrue = nTrue,
                                               printProgress = TRUE)

  return(simStudyResults)
}

set.seed(31415)
simStudyResults <- lapply(seq(from = 1, to = 16, by = 3), analysis_given_sidelength,
                          rmax = 3)

save(simStudyResults, file = "studyWindow.RDa")
