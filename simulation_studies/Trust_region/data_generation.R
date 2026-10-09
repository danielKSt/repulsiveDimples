source(file = "simulation_studies/Trust_region/Helper_scripts/study_helpers.R")

# With R's default generator set.seed does not make pmclapply reproducible, since the
# forked workers seed themselves. L'Ecuyer-CMRG gives each worker a stream derived from the
# seed instead. simulate_for_study_settings reseeds with the same seed on every call, so the
# thinning types are also thinned from the same Thomas patterns.
RNGkind("L'Ecuyer-CMRG")

studySettings <- data.frame(studyNr = c(1:3),
                             kappa = c(0.1234, 0.1234, 0.4321),
                             omega = c(0.8, 0.8, 0.4),
                             mu = c(2/0.1234, 2/0.1234, 4/0.4321),
                             rRange = c(0.4, 0.1, 0.2),
                             sidelength = c(40, 40, 30))
nTrue <- 1000000
nUnions <- 1000

for(i in studySettings$studyNr){
  setting <- studySettings[studySettings$studyNr == i, ]
  parThomas <- setting[, c("kappa", "omega", "mu", "rRange")]
  print(paste0("Simulating for study setting: ", i))

  print("Thinning type 1:")
  simulate_for_study_settings(parThomas = parThomas, studyNr = i,
                              nTrue = nTrue, nUnions = nUnions, sidelength = setting$sidelength,
                              thinningType = 1, mc.cores = 6)

  # commented out for thinningType = 2, since we have already simulated this
  #print("Thinning type 2:")
  # simulate_for_study_settings(parThomas = parThomas, studyNr = i,
  #                             nTrue = nTrue, nUnions = nUnions, sidelength = setting$sidelength,
  #                             thinningType = 2, mc.cores = 6)

  print("Thinning type 3:")
  simulate_for_study_settings(parThomas = parThomas, studyNr = i,
                              nTrue = nTrue, nUnions = nUnions, sidelength = setting$sidelength,
                              thinningType = 3, mc.cores = 6)
}
