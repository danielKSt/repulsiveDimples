# Targets for the Variance gamma trust region studies, the twin of
# simulation_studies/Thomas/Trust_region/data_generation.R.
#
# Studies 1-3 keep the Thomas studies' kappa, mu and repulsion range, with nu = 1/2 (the
# exponential kernel) and scale = omega/sqrt(3); studies 4-6 repeat them with nu = 5/2 and
# scale = omega/sqrt(7). A daughter's mean squared distance from its parent is
# 4*scale^2*(nu + 1), so both give 2*omega^2, as for the Thomas process, and the clusters
# are the same size as there. The nu = 1/2 mapping is also that of the Variance gamma
# importance sampling study 4. Every study uses a 30 x 30 window, and data are made for
# all three thinning types.
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
                            nu = c(1/2, 1/2, 1/2),
                            rRange = c(0.4, 0.1, 0.2),
                            sidelength = c(30, 30, 30))
# The nu = 5/2 studies need numbers of their own, since the data are saved by study number.
studySettings <- rbind(studySettings,
                       data.frame(studyNr = c(4:6),
                                  kappa = c(0.1234, 0.1234, 0.4321),
                                  scale = thomasOmega/sqrt(7),
                                  mu = c(2/0.1234, 2/0.1234, 4/0.4321),
                                  nu = c(5/2, 5/2, 5/2),
                                  rRange = c(0.4, 0.1, 0.2),
                                  sidelength = c(30, 30, 30)))
thinningTypes <- c(1:3)

if(dry_run){
  # The parent of tempdir() is the per-user temporary folder, which outlives this R
  # session, so the quad_run.R dry run can find what the data_generation.R one saved.
  studyRoot <- file.path(dirname(tempdir()), "varGamma_trust_region_dry_run/")
  nTrue <- 200
  nUnions <- 50
  nCores <- 4
} else {
  nTrue <- 1000000
  nUnions <- 1000
  # The patterns depend on the number of cores as well as on the seed, since every forked
  # worker draws from its own L'Ecuyer stream, so keep it fixed to reproduce the data.
  nCores <- 16
}
cat(sprintf("Variance gamma trust region data: studies %s, nTrue %d, nUnions %d, thinning types %s, %d cores, saving under %s%s\n",
            paste(studySettings$studyNr, collapse = ", "), nTrue, nUnions,
            paste(thinningTypes, collapse = ", "), nCores, studyRoot,
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
                                thinningType = thinningType, mc.cores = nCores)
  }
}
