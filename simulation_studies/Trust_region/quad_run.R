library(spatstat)
library(parallel)
library(mcprogress)
library(repulsiveDimples)

source("simulation_studies/Trust_region/Helper_scripts/study_run_helpers.R")

# Iterations allowed inside each pre-fitting block, which differs by study.
maxIterPrefit <- c(10, 10, 15)

run_quad <- function(studyNr, thinningType){
  # Loaded into an environment of its own, so nothing from one data file can be picked up
  # by the next study's fit.
  dat <- new.env()
  load(file = paste0(study_folder(thinningType), "Data/study", studyNr, ".RDa"), envir = dat)
  start_params <- study_start_params(K_hat_unions = dat$K_hat_unions,
                                     par_thomas = dat$par_thomas, rho_hat = dat$rho_hat)
  run_study(study = paste0("study", studyNr), variant = "quad", method = "quadratic",
            start_params = start_params, par_thomas = dat$par_thomas,
            rho_hat = dat$rho_hat, K_hat = dat$K_hat, thinningType = thinningType,
            max.iter.prefit = maxIterPrefit[studyNr])
  # run_study has already saved every level, so there is no need to keep its fits in memory.
  invisible(NULL)
}

studyNrs <- 1:3
# Type 2 is left out, since its results are already in Thinning_type2/Results and a rerun
# would overwrite them.
thinningTypes <- c(1, 3)

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
