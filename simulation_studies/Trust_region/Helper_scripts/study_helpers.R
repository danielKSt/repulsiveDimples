# Shared helpers for the trust region studies, in two sections: data generation, used by
# data_generation.R and combined_run_sample_error.R, and study runs, used by quad_run.R,
# combined_run_sample_error.R and the study1_conjugate.R and tolerance_pilot_run.R scripts
# in Thinning_type2/Study_scripts.
#
# The runs differ only in which data set they load, how many iterations the pre-fitting
# blocks get, which trust step method they ask for and which thinning type they fit, so
# everything else lives here rather than in copies that have to be kept in step by hand.

# Load libraries
suppressPackageStartupMessages({
  library(spatstat)
  library(ggplot2)
  library(parallel)
  library(mcprogress)
  library(repulsiveDimples)
})

# Each thinning type keeps its data and results in its own Thinning_type<N> folder. The
# file names are the same across types, so sharing one folder would let a run of one type
# silently overwrite another's.
study_folder <- function(thinningType){
  paste0("simulation_studies/Trust_region/Thinning_type", thinningType, "/")
}


# Data generation: ----

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

# The targets a fit is made against, rho_hat and K_hat from nTrue patterns, and K_hat_unions
# from nUnions further patterns, which only the starting parameters are derived from.
#
# With mc.cores = 1 the patterns are simulated serially with lapply, so the function can be
# called inside a worker that is itself one of many running in parallel, and the seed set
# in that worker governs the whole simulation.
estimate_targets <- function(parThomas, nTrue, nUnions, sidelength, thinningType = 2,
                             mc.cores = 6){
  map <- function(X, FUN, ...){
    if(mc.cores > 1){
      mcprogress::pmclapply(X = X, FUN = FUN, ..., mc.cores = mc.cores)
    } else {
      lapply(X = X, FUN = FUN, ...)
    }
  }

  sim_list_thomas <- map(X = rep(sidelength, nTrue), FUN = generate_data,
                         sim_pars = parThomas, r_vec = seq(from = 0, to = 3, by = 0.05),
                         thinningType = thinningType)

  rho_baseline <- sapply(sim_list_thomas, function(res) res$rho)
  rho_hat <- mean(rho_baseline)

  K_lambda_baseline <- lapply(sim_list_thomas, function(res) res$K)
  K_hat <- K_importance_sampling(w_is = rep(1, nTrue), K_lambda_baseline = K_lambda_baseline,
                                 rho_baseline = rho_baseline)

  # Want to also test with union approach:
  sim_list_thomas <- map(X = rep(parThomas$kappa, nUnions), FUN = rThomas_matern_thinned,
                         scale = parThomas$omega,
                         mu = parThomas$mu,
                         repulsionRange = parThomas$rRange,
                         xlims = c(0, sidelength),
                         ylims = c(0, sidelength),
                         saveparents = TRUE,
                         thinningType = thinningType)
  points_input <- lapply(X = sim_list_thomas, function(a) {
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

simulate_for_study_settings <- function(parThomas, studyNr, nTrue, nUnions, sidelength, thinningType = 2,
                                        mc.cores = 6, seed = 1350){
  set.seed(seed = seed)
  targets <- estimate_targets(parThomas = parThomas, nTrue = nTrue, nUnions = nUnions,
                              sidelength = sidelength, thinningType = thinningType,
                              mc.cores = mc.cores)
  rho_hat <- targets$rho_hat
  K_hat <- targets$K_hat
  K_hat_unions <- targets$K_hat_unions

  # Saved as par_thomas, the name run_study and the existing data files use.
  par_thomas <- parThomas
  save(par_thomas, rho_hat, K_hat, K_hat_unions, nTrue, nUnions, sidelength,
       file = paste0(study_folder(thinningType), "Data/study", studyNr, ".RDa"))
}


# Study runs: ----

# Starting parameters: fit as if the pattern were unthinned, and back out mu from the
# intensity the thinning would have produced.
# We use the intensity for thinningType 2 for all cases to avoid Lambert's W function in type 1,
# and since type 3 is not available analyticaly
study_start_params <- function(K_hat_unions, par_thomas, rho_hat){
  start_fitted <- thomas.estK(X = K_hat_unions, rmin = 2*par_thomas$rRange, rmax = 3)$par
  lambda <- -log(1 - pi*rho_hat*par_thomas$rRange^2)/(pi*par_thomas$rRange^2)
  start_params <- log(c(start_fitted[1], start_fitted[2], lambda/start_fitted[1]))
  names(start_params) <- c("kappa", "omega", "mu")
  return(start_params)
}

# One timed fit, with the settings every trust region study shares.
#
# method         "conjugate" or "quadratic", passed to both phases of the fit
# max.iter.prefit iterations allowed inside each pre-fitting block, which differs by study
# tolPrefitLoops NULL leaves it out of the fit call, so the package default applies
#
# The seed is not set here: the caller sets it, and should also have set
# options(mc.cores = 1), see run_study for why.
fit_trust_region <- function(start_params, par_thomas, rho_hat, K_hat, nSims, thinningType,
                             method, max.iter.prefit, tolPrefitLoops = NULL){
  common <- list(params = start_params, parFreeIndex = c(1, 2, 3),
                 repRange = par_thomas$rRange, rho_hat = rho_hat, K_hat = K_hat,
                 xlims = c(0, 3), ylims = c(0, 3), nSims = nSims, thinningType = thinningType,
                 deltaInit = 0.2, eta = 0.05, deltaMax = 1.0, deltaMin = 0.0001,
                 wq = c(1000, 1/4), normalized = TRUE, eta_trust = 0.6,
                 eta_converged = 0.999, validate.converged = FALSE,
                 maxPrefitLoops = 8, max.iter.prefit = max.iter.prefit,
                 tol = 10^-8, max.iter = 30, printProgress = FALSE)
  # NULL leaves tolPrefitLoops out of the call entirely, so the package default applies
  # and the study scripts keep saying "default" rather than pinning a number that would
  # then have to be chased if the default ever moved. The tolerance pilot passes a value.
  if(!is.null(tolPrefitLoops)){
    common$tolPrefitLoops <- tolPrefitLoops
  }
  # Only the arguments the chosen method actually reads are passed, so a script does not
  # imply it tuned something the method ignores. The quadratic step's own two knobs,
  # interp_fraction and bisection_iterations, are left at their defaults.
  perMethod <- if(method == "conjugate"){
    list(method.main = "conjugate", method.prefit = "conjugate",
         carry.directions = TRUE,
         subsection_count.prefit = 8, line_iterations.prefit = 3,
         max.directions.prefit = NULL,
         subsection_count.main = 11, line_iterations.main = 4,
         max.directions.main = NULL)
  } else {
    list(method.main = "quadratic", method.prefit = "quadratic")
  }
  # Timed inside the worker and around the fit alone, so neither the fork nor the
  # scheduling is counted. Both clocks are kept because they answer different questions.
  # The fits run nCores at a time, so `secs` is what the study costs in practice but is
  # inflated by whatever the sibling workers are doing; `cpu_secs` is the work this fit
  # actually did, and is the fairer number for comparing the two trust step methods.
  # They coincide only when nCores is well under the machine's core count.
  #
  # The times go into the returned list rather than around it, so everything downstream
  # that reads a fit by name -- getParamsFinal and the rest of
  # results_analysis_helpers.R -- is unaffected.
  t0  <- proc.time()
  res <- do.call(min_contrast_trust_region, c(common, perMethod))
  el  <- proc.time() - t0
  res$secs     <- unname(el[["elapsed"]])
  res$cpu_secs <- unname(el[["user.self"]] + el[["sys.self"]])
  return(res)
}

# One study at one ensemble size, over nRuns fits.
#
# study          "study1", "study2", "study3" -- names the data and the output files
# method         "conjugate" or "quadratic", passed to both phases of the fit
# variant        the tag the output files carry, "conjugate" or "quad"
# max.iter.prefit iterations allowed inside each pre-fitting block, which differs by study
# thinningType   the Matern thinning type, 1, 2 or 3. It also decides the folder the results
#                are saved in, see study_folder, so it has no default: a script that left
#                it out would fit the wrong type and file the results under another type.
#
# Every fit is timed individually: `secs` and `cpu_secs` are added to the fit object it
# returns, so a saved level holds one pair of times per fit rather than only the level
# total.
#
# tolPrefitLoops is NULL by default, which leaves it out of the fit call so the package
# default applies. Pass a number to override it, as tolerance_pilot_run.R does.
#
# Pre-fitting stops on tolPrefitLoops, left at the package default, with maxPrefitLoops as
# a cap that should rarely bind. That is the point of running these: a fixed number of
# cycles cannot know where the blocks stop making progress, and the two studies differ in
# how quickly they get there.
run_study <- function(study, variant, method, start_params, par_thomas, rho_hat, K_hat,
                      max.iter.prefit, thinningType, nRuns = 200, nCores = 32,
                      nSimsLevels = c(500, 1000, 5000, 10000),
                      seedBases = c(10000, 20000, 30000, 40000),
                      tolPrefitLoops = NULL){
  fit_one <- function(i, nSims, seedBase){
    # Inner loops serial, at every nSims, for two reasons.
    #
    # Speed: parallel importance sampling is a net loss below nSims = 1000 -- measured at
    # 0.57x for nSims = 500, fork overhead beating the per-pattern work on patterns this
    # small -- and even at nSims = 10000 it only reaches 1.66x. Those figures came off a
    # 10-core machine, but the conclusion does not depend on the core count: inner
    # parallelism buys at most about 1.7x however many cores it takes, while the same
    # cores running whole fits alongside each other scale roughly linearly. Spend them on
    # the outer loop instead.
    #
    # Reproducibility: with one core mclapply never forks, so the seed set below governs
    # the whole fit. Leave it out and the inner forks seed themselves, and the run cannot
    # be reproduced. The seeds match across the two methods, so a conjugate fit and a
    # quadratic fit with the same index see the same ensembles and are directly paired.
    options(mc.cores = 1)
    set.seed(seedBase + i)
    fit_trust_region(start_params = start_params, par_thomas = par_thomas,
                     rho_hat = rho_hat, K_hat = K_hat, nSims = nSims,
                     thinningType = thinningType, method = method,
                     max.iter.prefit = max.iter.prefit, tolPrefitLoops = tolPrefitLoops)
  }

  # Each level is saved as it finishes. nSims = 10000 is hours of work on its own, and
  # mclapply hands back a "try-error" for a worker that died rather than aborting the
  # rest, so a single bad fit costs one fit instead of the study. The level files are only
  # a safeguard while the study runs: once every level is done and the combined file is
  # saved, they are deleted.
  run_level <- function(nSims, seedBase){
    print(paste0("Starting ", study, " ", variant, " run with nSims = ", nSims))
    t0  <- Sys.time()
    res <- mcprogress::pmclapply(seq_len(nRuns), fit_one, nSims = nSims, seedBase = seedBase,
                                 mc.cores = nCores, mc.preschedule = FALSE)
    save(res, nRuns, nSims, seedBase, method, thinningType, start_params, par_thomas,
         rho_hat, K_hat, file = level_file(nSims))
    ok <- !vapply(res, inherits, logical(1), "try-error")
    # A level in which every fit died has no times to summarise, and min/max of nothing
    # would turn the report into +/-Inf. Report what there is.
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
  # One element per ensemble size, each list(res = <fits>, nSims = <size>), which is the
  # shape combineParamsFinal in results_analysis_helpers.R expects. seedBases is saved
  # since the level files, which also record each level's seed, are deleted below.
  save(res, nRuns, method, thinningType, seedBases, start_params, par_thomas, rho_hat, K_hat,
       file = paste0(resFolder, study, "_", variant, ".RDa"))
  # save() stops with an error if the combined file cannot be written, so the level files
  # are only removed once everything in them is in the combined file.
  file.remove(level_file(unique(nSimsLevels)))
  return(invisible(res))
}

# One replicate of the sample error study, combined_run_sample_error.R: a fresh sample of
# nTrue and nUnions patterns, the targets and starting parameters estimated from it, and one
# fit per value of nSims, all on the same sample.
#
# i             the replicate's index, added to the seed bases
# setting       one row of the study settings, with the Thomas parameters, rRange,
#               sidelength, nTrue and nUnions
# maxIterPrefit iterations allowed inside each pre-fitting block
# nSims         the ensemble sizes to fit at, one fit each
# dataSeedBase  the sample is simulated from dataSeedBase + i
# fitSeedBases  one per value of nSims, the fit at nSims[j] uses fitSeedBases[j] + i
#
# Simulation and fits are serial inside the worker, so the seeds set here govern all of it,
# see run_study.
single_sample <- function(i, setting, thinningType, maxIterPrefit, nSims, dataSeedBase,
                          fitSeedBases){
  options(mc.cores = 1)
  parThomas <- setting[, c("kappa", "omega", "mu", "rRange")]

  set.seed(dataSeedBase + i)
  sample <- try({
    targets <- estimate_targets(parThomas = parThomas, nTrue = setting$nTrue,
                                nUnions = setting$nUnions, sidelength = setting$sidelength,
                                thinningType = thinningType, mc.cores = 1)
    # thomas.estK can fail, or give an unusable start, on K estimated from this few patterns.
    startParams <- study_start_params(K_hat_unions = targets$K_hat_unions,
                                      par_thomas = parThomas, rho_hat = targets$rho_hat)
    if(any(!is.finite(startParams))){
      stop("non-finite starting parameters")
    }
    list(rho_hat = targets$rho_hat, K_hat = targets$K_hat, start_params = startParams)
  }, silent = TRUE)
  # A failed sample has nothing to fit, so each of its fits is recorded as the error, which
  # combine_params_final drops like any other failed fit.
  if(inherits(sample, "try-error")){
    return(list(sample = sample, fits = rep(list(sample), length(nSims))))
  }

  fits <- lapply(seq_along(nSims), function(j){
    set.seed(fitSeedBases[j] + i)
    try(fit_trust_region(start_params = sample$start_params, par_thomas = parThomas,
                         rho_hat = sample$rho_hat, K_hat = sample$K_hat, nSims = nSims[j],
                         thinningType = thinningType, method = "quadratic",
                         max.iter.prefit = maxIterPrefit), silent = TRUE)
  })
  return(list(sample = sample, fits = fits))
}
