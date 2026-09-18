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
# Pre-fitting stops on tolPrefitLoops, left at the package default, with maxPrefitLoops as
# a cap that should rarely bind. That is the point of running these: a fixed number of
# cycles cannot know where the blocks stop making progress, and the two studies differ in
# how quickly they get there.
run_study <- function(study, variant, method, start_params, par_thomas, rho_hat, K_hat,
                      max.iter.prefit, nRuns = 200, nCores = 32,
                      nSimsLevels = c(500, 1000, 5000, 10000),
                      seedBases = c(10000, 20000, 30000, 40000)){
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
                   # tolPrefitLoops left at the package default on purpose.
                   maxPrefitLoops = 8, max.iter.prefit = max.iter.prefit,
                   tol = 10^-8, max.iter = 30, printProgress = FALSE)
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
    do.call(min_contrast_trust_region, c(common, perMethod))
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
    message(sprintf("%s %s, nSims = %5d: %.1f min, %d of %d fits failed",
                    study, variant, nSims,
                    as.numeric(difftime(Sys.time(), t0, units = "mins")),
                    sum(vapply(res, inherits, logical(1), "try-error")), nRuns))
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
