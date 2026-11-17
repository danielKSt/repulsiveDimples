# Shared helpers for the Variance gamma trust region studies, the twin of
# simulation_studies/Thomas/Trust_region/Helper_scripts/study_helpers.R, in two sections:
# data generation, used by data_generation.R, and study runs, used by quad_run.R.
#
# Kept as a separate copy rather than sourced from the Thomas helpers, as for the
# importance sampling studies: the simulator, the starting values and the fit function all
# differ, and the Thomas copy is the provenance of results that are already saved.
#
# Parameters are (kappa, scale, mu, nu), with nu fixed at the value the data were
# simulated with. There is no pre-fitting, see min_contrast_trust_region_varGamma.

# Load libraries
suppressPackageStartupMessages({
  library(spatstat)
  library(ggplot2)
  library(parallel)
  library(mcprogress)
  library(repulsiveDimples)
})

# Each thinning type keeps its data and results in its own Thinning_type<N> folder, as for
# the Thomas studies. studyRoot can be pointed elsewhere after sourcing, which the
# --dry-run of the run scripts does so that a test run never overwrites real data.
studyRoot <- "simulation_studies/VarGamma/Trust_region/"
study_folder <- function(thinningType){
  paste0(studyRoot, "Thinning_type", thinningType, "/")
}


# Data generation: ----

generate_data <- function(sidelength, sim_pars, r_vec, thinningType = 2){
  pattern <- rVarGamma_matern_thinned(kappa = sim_pars$kappa,
                                      scale = sim_pars$scale,
                                      mu = sim_pars$mu,
                                      nu = sim_pars$nu,
                                      repulsionRange = sim_pars$rRange,
                                      xlims = c(0, sidelength),
                                      ylims = c(0, sidelength),
                                      thinningType = thinningType,
                                      saveparents = TRUE)
  rho_est <- estimate_rho_baseline(pattern = pattern)
  K_est <- estimate_K_lambda_baseline(pattern = pattern, r_vec = r_vec)
  return(list(rho = rho_est, K = K_est))
}

# The targets a fit is made against, rho_hat and K_hat from nTrue patterns, and K_hat_unions
# from nUnions further patterns, which only the starting parameters are derived from.
#
# With mc.cores = 1 the patterns are simulated serially with lapply, so the function can be
# called inside a worker that is itself one of many running in parallel, and the seed set
# in that worker governs the whole simulation.
estimate_targets <- function(parVarGamma, nTrue, nUnions, sidelength, thinningType = 2,
                             mc.cores = 6){
  map <- function(X, FUN, ...){
    if(mc.cores > 1){
      mcprogress::pmclapply(X = X, FUN = FUN, ..., mc.cores = mc.cores)
    } else {
      lapply(X = X, FUN = FUN, ...)
    }
  }

  sim_list <- map(X = rep(sidelength, nTrue), FUN = generate_data,
                  sim_pars = parVarGamma, r_vec = seq(from = 0, to = 3, by = 0.05),
                  thinningType = thinningType)

  rho_baseline <- sapply(sim_list, function(res) res$rho)
  rho_hat <- mean(rho_baseline)

  K_lambda_baseline <- lapply(sim_list, function(res) res$K)
  K_hat <- K_importance_sampling(w_is = rep(1, nTrue), K_lambda_baseline = K_lambda_baseline,
                                 rho_baseline = rho_baseline)

  # The union approach, for the starting values:
  sim_list <- map(X = rep(parVarGamma$kappa, nUnions), FUN = rVarGamma_matern_thinned,
                  scale = parVarGamma$scale,
                  mu = parVarGamma$mu,
                  nu = parVarGamma$nu,
                  repulsionRange = parVarGamma$rRange,
                  xlims = c(0, sidelength),
                  ylims = c(0, sidelength),
                  saveparents = TRUE,
                  thinningType = thinningType)
  points_input <- lapply(X = sim_list, function(a) {
    if(nrow(a$thinned) == 0){
      return(0)
    } else {
      return(a$thinned)
    }
  })
  K_hat_unions <- K_est.unions(points_input = points_input, l = sidelength, spacing = 5,
                               r_vec = seq(from = 0, to = 3, by = 0.05), timescale = 1)

  return(list(rho_hat = rho_hat, K_hat = K_hat, K_hat_unions = K_hat_unions))
}

simulate_for_study_settings <- function(parVarGamma, studyNr, nTrue, nUnions, sidelength,
                                        thinningType = 2, mc.cores = 6, seed = 1350){
  set.seed(seed = seed)
  targets <- estimate_targets(parVarGamma = parVarGamma, nTrue = nTrue, nUnions = nUnions,
                              sidelength = sidelength, thinningType = thinningType,
                              mc.cores = mc.cores)
  rho_hat <- targets$rho_hat
  K_hat <- targets$K_hat
  K_hat_unions <- targets$K_hat_unions

  par_varGamma <- parVarGamma
  dataFolder <- paste0(study_folder(thinningType), "Data/")
  dir.create(dataFolder, showWarnings = FALSE, recursive = TRUE)
  save(par_varGamma, rho_hat, K_hat, K_hat_unions, nTrue, nUnions, sidelength,
       file = paste0(dataFolder, "study", studyNr, ".RDa"))
}


# Study runs: ----

# Starting parameters: fit the unthinned Variance gamma K-function to the union estimate
# with nu fixed at its true value, and back out mu from the intensity the thinning would
# have produced. As for the Thomas studies, the type 2 intensity is used for every thinning
# type, to avoid Lambert's W function in type 1, and since type 3 has no closed form.
#
# Returns (log kappa, log scale, log mu, nu), the parameter vector
# min_contrast_trust_region_varGamma takes.
study_start_params <- function(K_hat_unions, par_varGamma, rho_hat){
  # $clustpar rather than $par, which calls the scale eta (Jalilian et al.'s name for it).
  start_fitted <- vargamma.estK(X = K_hat_unions, nu = par_varGamma$nu,
                                rmin = 2*par_varGamma$rRange, rmax = 3)$clustpar
  lambda <- -log(1 - pi*rho_hat*par_varGamma$rRange^2)/(pi*par_varGamma$rRange^2)
  start_params <- c(log(c(start_fitted[["kappa"]], start_fitted[["scale"]],
                          lambda/start_fitted[["kappa"]])),
                    par_varGamma$nu)
  names(start_params) <- c("kappa", "scale", "mu", "nu")
  return(start_params)
}

# One timed fit, with the settings every Variance gamma trust region study shares. These
# are the Thomas studies' main fit settings, except that eta_trust is an argument.
#
# The seed is not set here: the caller sets it, and should also have set
# options(mc.cores = 1), see run_study for why.
fit_trust_region <- function(start_params, par_varGamma, rho_hat, K_hat, nSims, thinningType,
                             eta_trust){
  # Timed inside the worker and around the fit alone. `secs` is what the study costs in
  # practice, inflated by whatever the sibling workers are doing, and `cpu_secs` is the
  # work this fit actually did.
  t0  <- proc.time()
  res <- min_contrast_trust_region_varGamma(params = start_params, parFreeIndex = c(1, 2, 3),
                                            repRange = par_varGamma$rRange,
                                            rho_hat = rho_hat, K_hat = K_hat,
                                            xlims = c(0, 3), ylims = c(0, 3), nSims = nSims,
                                            thinningType = thinningType,
                                            deltaInit = 0.2, eta = 0.05, deltaMax = 1.0,
                                            deltaMin = 0.0001, wq = c(1000, 1/4),
                                            normalized = TRUE, eta_trust = eta_trust,
                                            eta_converged = 0.999, validate.converged = FALSE,
                                            method = "quadratic",
                                            tol = 10^-8, max.iter = 30, printProgress = FALSE)
  el  <- proc.time() - t0
  res$secs     <- unname(el[["elapsed"]])
  res$cpu_secs <- unname(el[["user.self"]] + el[["sys.self"]])
  return(res)
}

# One study at every ensemble size in nSimsLevels, over nRuns fits each.
#
# study          "study1", "study2", "study3" -- names the data and the output files
# variant        the tag the output files carry
# thinningType   the Matern thinning type, 1, 2 or 3. It also decides the folder the results
#                are saved in, so it has no default.
# eta_trust      passed to the fit, see min_contrast_trust_region_varGamma
#
# Each level is saved as it finishes, and the level files are deleted once the combined
# file is saved, as for the Thomas studies.
run_study <- function(study, variant, start_params, par_varGamma, rho_hat, K_hat,
                      thinningType, eta_trust, nRuns = 200, nCores = 32,
                      nSimsLevels = c(500, 1000, 5000, 10000),
                      seedBases = c(10000, 20000, 30000, 40000)){
  fit_one <- function(i, nSims, seedBase){
    # Inner loops serial, so the cores go to whole fits running alongside each other, and
    # so the seed set below governs the whole fit. See the Thomas run_study for the
    # measurements behind this.
    options(mc.cores = 1)
    set.seed(seedBase + i)
    fit_trust_region(start_params = start_params, par_varGamma = par_varGamma,
                     rho_hat = rho_hat, K_hat = K_hat, nSims = nSims,
                     thinningType = thinningType, eta_trust = eta_trust)
  }

  run_level <- function(nSims, seedBase){
    print(paste0("Starting ", study, " ", variant, " run with nSims = ", nSims))
    t0  <- Sys.time()
    res <- mcprogress::pmclapply(seq_len(nRuns), fit_one, nSims = nSims, seedBase = seedBase,
                                 mc.cores = nCores, mc.preschedule = FALSE)
    save(res, nRuns, nSims, seedBase, thinningType, eta_trust, start_params, par_varGamma,
         rho_hat, K_hat, file = level_file(nSims))
    ok <- !vapply(res, inherits, logical(1), "try-error")
    timing <- if(any(ok)){
      perFit <- vapply(res[ok], function(f) f$secs, numeric(1))
      perCpu <- vapply(res[ok], function(f) f$cpu_secs, numeric(1))
      sprintf("; per fit median %.1f s (cpu %.1f s), range %.1f-%.1f s",
              stats::median(perFit), stats::median(perCpu), min(perFit), max(perFit))
    } else {
      ""
    }
    message(sprintf("%s %s, nSims = %5d: %.1f min, %d of %d fits failed%s",
                    study, variant, nSims,
                    as.numeric(difftime(Sys.time(), t0, units = "mins")),
                    sum(!ok), nRuns, timing))
    return(list(res = res, nSims = nSims))
  }

  resFolder <- paste0(study_folder(thinningType), "Results/")
  dir.create(resFolder, showWarnings = FALSE, recursive = TRUE)
  level_file <- function(nSims){
    paste0(resFolder, study, "_", variant, "_nSims", nSims, ".RDa")
  }

  res <- vector(mode = "list", length = length(nSimsLevels))
  for(i in seq_along(nSimsLevels)){
    res[[i]] <- run_level(nSims = nSimsLevels[i], seedBase = seedBases[i])
  }
  # One element per ensemble size, each list(res = <fits>, nSims = <size>), the shape the
  # Thomas results_analysis_helpers.R expects.
  save(res, nRuns, thinningType, eta_trust, seedBases, start_params, par_varGamma, rho_hat,
       K_hat, file = paste0(resFolder, study, "_", variant, ".RDa"))
  file.remove(level_file(unique(nSimsLevels)))
  return(invisible(res))
}
