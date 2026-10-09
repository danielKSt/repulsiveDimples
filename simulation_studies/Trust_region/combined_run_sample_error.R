source("simulation_studies/Trust_region/Helper_scripts/study_helpers.R")
RNGkind("L'Ecuyer-CMRG")


studySettings <- data.frame(studyNr = c(1:3),
                            kappa = c(0.1234, 0.1234, 0.4321),
                            omega = c(0.8, 0.8, 0.4),
                            mu = c(2/0.1234, 2/0.1234, 4/0.4321),
                            rRange = c(0.4, 0.1, 0.2),
                            sidelength = c(20, 15, 15),
                            nTrue = c(10, 10, 10),
                            nUnions = c(10, 10, 10))

# Iterations allowed inside each pre-fitting block, which differs by study, as in quad_run.R.
maxIterPrefit <- c(10, 10, 15)

nResims <- 200
nSims <- c(500, 1000, 2000)
nCores <- 32
thinningTypes <- 1:3

# Seeds are the same across studies and thinning types, as in run_study. The data of
# replicate i is simulated from dataSeedBase + i, and its fit at nSims[j] from
# fitSeedBases[j] + i. With nResims <= 9999 none of them overlap.
dataSeedBase <- 90000
fitSeedBases <- c(10000, 20000, 30000)

for(thinningType in thinningTypes){
  print(paste0("Thinning type: ", thinningType))
  dir.create(paste0(study_folder(thinningType), "Results/"), showWarnings = FALSE,
             recursive = TRUE)
  for(i in studySettings$studyNr){
    resFile <- paste0(study_folder(thinningType), "Results/study", i, "_sample_error.RDa")
    # A finished combination is not run again, so the script can be restarted after a
    # crash. Delete the file to rerun it.
    if(file.exists(resFile)){
      message("Skipping study ", i, " with thinning type ", thinningType, ", ", resFile,
              " already exists")
      next
    }
    print(paste0("Study setting: ", i))
    setting <- studySettings[studySettings$studyNr == i, ]
    par_thomas <- setting[, c("kappa", "omega", "mu", "rRange")]

    t0 <- Sys.time()
    reps <- mcprogress::pmclapply(seq_len(nResims), single_sample, setting = setting,
                                  thinningType = thinningType,
                                  maxIterPrefit = maxIterPrefit[i],
                                  nSims = nSims, dataSeedBase = dataSeedBase,
                                  fitSeedBases = fitSeedBases,
                                  mc.cores = nCores, mc.preschedule = FALSE)

    # A worker that died returns a try-error rather than the replicate's list.
    reps <- lapply(reps, function(r){
      if(inherits(r, "try-error")) list(sample = r, fits = rep(list(r), length(nSims))) else r
    })
    # Rearranged to one element per ensemble size, each list(res = <fits>, nSims = <size>),
    # the shape combine_params_final in results_analysis_helpers.R expects. Fit k at every
    # level was made on sample k, so the levels are paired by index.
    res <- lapply(seq_along(nSims), function(j){
      list(res = lapply(reps, function(r) r$fits[[j]]), nSims = nSims[j])
    })
    samples <- lapply(reps, function(r) r$sample)

    nFailedSamples <- sum(vapply(samples, inherits, logical(1), "try-error"))
    nFailedFits <- sum(vapply(res, function(level)
      sum(vapply(level$res, inherits, logical(1), "try-error")), numeric(1)))
    message(sprintf("Study %d, thinning type %d: %.1f min, %d of %d samples failed, %d of %d fits failed",
                    i, thinningType, as.numeric(difftime(Sys.time(), t0, units = "mins")),
                    nFailedSamples, nResims, nFailedFits, nResims*length(nSims)))

    save(res, samples, setting, par_thomas, thinningType, nResims, nSims,
         dataSeedBase, fitSeedBases, file = resFile)
    rm(reps, res, samples)
  }
}
