# Trust region fits for the Variance gamma studies, the twin of
# simulation_studies/Thomas/Trust_region/quad_run.R: 200 fits per study at each of
# nSims = 500, 1000, 5000 and 10000, with the quadratic trust step and no pre-fitting.
# Needs the data from data_generation.R.
#
# Run from the project root. The results are saved in
# simulation_studies/VarGamma/Trust_region/Thinning_type<N>/Results.
#
# For a fast end-to-end check, add --dry-run. It reads the data the data_generation.R dry
# run saved in the temporary folder, so run that first:
#   Rscript simulation_studies/VarGamma/Trust_region/data_generation.R --dry-run
#   Rscript simulation_studies/VarGamma/Trust_region/quad_run.R --dry-run

source("simulation_studies/VarGamma/Trust_region/Helper_scripts/study_helpers.R")

dry_run <- "--dry-run" %in% commandArgs(trailingOnly = TRUE)

# Smallest effective sample size, as a fraction of nSims, a proposed iterate may have.
# The Thomas studies used 0.6.
eta_trust <- 0.6

studyNrs <- 1:3
thinningTypes <- 2

if(dry_run){
  # The parent of tempdir() is the per-user temporary folder, which outlives this R
  # session, so the quad_run.R dry run can find what the data_generation.R one saved.
  studyRoot <- file.path(dirname(tempdir()), "varGamma_trust_region_dry_run/")
  runSettings <- list(nRuns = 2, nCores = 2, nSimsLevels = c(50, 100), seedBases = c(10000, 20000))
} else {
  runSettings <- list(nRuns = 200, nCores = 32, nSimsLevels = c(500, 1000, 5000, 10000),
                      seedBases = c(10000, 20000, 30000, 40000))
}

run_quad <- function(studyNr, thinningType){
  # Loaded into an environment of its own, so nothing from one data file can be picked up
  # by the next study's fit.
  dat <- new.env()
  load(file = paste0(study_folder(thinningType), "Data/study", studyNr, ".RDa"), envir = dat)
  start_params <- study_start_params(K_hat_unions = dat$K_hat_unions,
                                     par_varGamma = dat$par_varGamma, rho_hat = dat$rho_hat)
  print(paste0("Starting parameters: ", paste(names(start_params), "=",
                                              signif(c(exp(start_params[1:3]), start_params[4]), 4),
                                              collapse = ", ")))
  do.call(run_study, c(list(study = paste0("study", studyNr), variant = "quad",
                            start_params = start_params, par_varGamma = dat$par_varGamma,
                            rho_hat = dat$rho_hat, K_hat = dat$K_hat,
                            thinningType = thinningType, eta_trust = eta_trust),
                       runSettings))
  # run_study has already saved every level, so there is no need to keep its fits in memory.
  invisible(NULL)
}

# Check that every data file is in place before starting, rather than finding a missing one
# a day into the run.
dataFiles <- outer(studyNrs, thinningTypes, function(studyNr, thinningType)
  paste0(study_folder(thinningType), "Data/study", studyNr, ".RDa"))
if(!all(file.exists(dataFiles))){
  stop("Missing data files: ", paste(dataFiles[!file.exists(dataFiles)], collapse = ", "))
}

for(studyNr in studyNrs){
  for(thinningType in thinningTypes){
    print(paste0("Fitting study setting ", studyNr, " with thinning type ", thinningType))
    run_quad(studyNr = studyNr, thinningType = thinningType)
  }
}
