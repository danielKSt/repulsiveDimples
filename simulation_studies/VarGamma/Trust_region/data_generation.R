# Targets for the Variance gamma trust region studies, the twin of
# simulation_studies/Thomas/Trust_region/data_generation.R.
#
# The three study settings keep the Thomas studies' kappa, mu, repulsion range and window,
# with nu fixed at 1/2 (the exponential kernel) and scale = omega/sqrt(3), so that the
# clusters are the same size as there: at nu = 1/2 a daughter's mean squared distance
# from its parent is 6*scale^2, against 2*omega^2 for the Thomas process. This is the same
# mapping as the Variance gamma importance sampling study 4.
#
# Run from the project root. The data are saved in
# simulation_studies/VarGamma/Trust_region/Thinning_type<N>/Data.
#
# For a fast end-to-end check, which saves to a temporary folder instead, add --dry-run:
#   Rscript simulation_studies/VarGamma/Trust_region/data_generation.R --dry-run

source(file = "simulation_studies/VarGamma/Trust_region/Helper_scripts/study_helpers.R")

dry_run <- "--dry-run" %in% commandArgs(trailingOnly = TRUE)

# With R's default generator set.seed does not make pmclapply reproducible, since the
# forked workers seed themselves. L'Ecuyer-CMRG gives each worker a stream derived from the
# seed instead. simulate_for_study_settings reseeds with the same seed on every call, so the
# thinning types are also thinned from the same Variance gamma patterns.
RNGkind("L'Ecuyer-CMRG")

thomasOmega <- c(0.8, 0.8, 0.4)
studySettings <- data.frame(studyNr = c(1:3),
                            kappa = c(0.1234, 0.1234, 0.4321),
                            scale = thomasOmega/sqrt(3),
                            mu = c(2/0.1234, 2/0.1234, 4/0.4321),
                            nu = 1/2,
                            rRange = c(0.4, 0.1, 0.2),
                            sidelength = c(40, 40, 30))
thinningTypes <- 2

if(dry_run){
  # The parent of tempdir() is the per-user temporary folder, which outlives this R
  # session, so the quad_run.R dry run can find what the data_generation.R one saved.
  studyRoot <- file.path(dirname(tempdir()), "varGamma_trust_region_dry_run/")
  nTrue <- 200
  nUnions <- 50
} else {
  nTrue <- 1000000
  nUnions <- 1000
}
cat(sprintf("Variance gamma trust region data: nTrue %d, nUnions %d, thinning types %s, saving under %s%s\n",
            nTrue, nUnions, paste(thinningTypes, collapse = ", "), studyRoot,
            if(dry_run) "   [DRY RUN]" else ""))

for(i in studySettings$studyNr){
  setting <- studySettings[studySettings$studyNr == i, ]
  parVarGamma <- setting[, c("kappa", "scale", "mu", "nu", "rRange")]
  print(paste0("Simulating for study setting: ", i))

  for(thinningType in thinningTypes){
    print(paste0("Thinning type ", thinningType, ":"))
    simulate_for_study_settings(parVarGamma = parVarGamma, studyNr = i,
                                nTrue = nTrue, nUnions = nUnions,
                                sidelength = setting$sidelength,
                                thinningType = thinningType, mc.cores = 6)
  }
}
