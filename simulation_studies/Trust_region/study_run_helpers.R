# Shared body of the study<n>_conjugate.R and study<n>_quad.R run scripts.
#
# Sourced by each of them. The six scripts differ only in which data set they load, how
# many iterations the pre-fitting blocks get, and which trust step method they ask for, so
# everything else lives here rather than in six copies that have to be kept in step by
# hand.

resFolder <- "simulation_studies/Trust_region/Results/"

# Starting parameters: fit as if the pattern were unthinned, and back out mu from the
# intensity the thinning would have produced.
study_start_params <- function(K_hat_unions, par_thomas, rho_hat){
  start_fitted <- thomas.estK(X = K_hat_unions, rmin = 2*par_thomas$rRange, rmax = 3)$par
  lambda <- -log(1 - pi*rho_hat*par_thomas$rRange^2)/(pi*par_thomas$rRange^2)
  start_params <- log(c(start_fitted[1], start_fitted[2], lambda/start_fitted[1]))
  names(start_params) <- c("kappa", "omega", "mu")
  return(start_params)
}

# One study at one ensemble size, over nRuns fits.
#
# study          "study1", "study2", "study3" -- names the data and the output files
# method         "conjugate" or "quadratic", passed to both phases of the fit
# variant        the tag the output files carry, "conjugate" or "quad"
# max.iter.prefit iterations allowed inside each pre-fitting block, which differs by study
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
                      max.iter.prefit, nRuns = 200, nCores = 32,
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

    common <- list(params = start_params, parFreeIndex = c(1, 2, 3),
                   repRange = par_thomas$rRange, rho_hat = rho_hat, K_hat = K_hat,
                   xlims = c(0, 3), ylims = c(0, 3), nSims = nSims,
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

  # Each level is saved as it finishes. nSims = 10000 is hours of work on its own, and
  # mclapply hands back a "try-error" for a worker that died rather than aborting the
  # rest, so a single bad fit costs one fit instead of the study.
  run_level <- function(nSims, seedBase){
    print(paste0("Starting ", study, " ", variant, " run with nSims = ", nSims))
    t0  <- Sys.time()
    res <- mcprogress::pmclapply(seq_len(nRuns), fit_one, nSims = nSims, seedBase = seedBase,
                                 mc.cores = nCores, mc.preschedule = FALSE)
    save(res, nRuns, nSims, seedBase, method, start_params, par_thomas, rho_hat, K_hat,
         file = paste0(resFolder, study, "_", variant, "_nSims", nSims, ".RDa"))
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

  res <- vector(mode = "list", length = length(nSimsLevels))
  for(i in seq_along(nSimsLevels)){
    res[[i]] <- run_level(nSims = nSimsLevels[i], seedBase = seedBases[i])
  }
  # One element per ensemble size, each list(res = <fits>, nSims = <size>), which is the
  # shape combineParamsFinal in results_analysis_helpers.R expects.
  save(res, nRuns, method, start_params, par_thomas, rho_hat, K_hat,
       file = paste0(resFolder, study, "_", variant, ".RDa"))
  return(invisible(res))
}
