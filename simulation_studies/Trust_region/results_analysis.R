suppressPackageStartupMessages({
  library(spatstat)
  library(ggplot2)
  library(repulsiveDimples)
  library(dplyr)
})

source(file = "simulation_studies/Trust_region/Helper_scripts/results_analysis_helpers.R")

thinningType <- 2
studyFolder <- paste0("simulation_studies/Trust_region/Thinning_type", thinningType, "/")

type <- c("quad", "sample_error")

type <- "sample_error"
type <- "quad"
studyNr <- 3

# Fetch start values for the parameters: ----
load(file = paste0(studyFolder, "Data/study", studyNr, ".RDa"))

startFitted <- thomas.estK(X = K_hat_unions, rmin = 2*par_thomas$rRange, rmax = 3)
lambda <- -log(1-pi*rho_hat*par_thomas$rRange^2)/(pi*par_thomas$rRange^2)
startFitted <- startFitted$par

startParams <- log(c(startFitted[1], startFitted[2], lambda/startFitted[1]))
names(startParams) <- c("kappa", "omega", "mu")
rm(startFitted, lambda, K_hat_unions, nTrue, sidelength, nUnions)
# The sample error fits each start from their own sample, so no single starting value applies.
if(type == "sample_error"){
  startParams <- NULL
}

# Look at the results: ----
load(file = paste0(studyFolder, "Results/study", studyNr, "_", type, ".RDa"))

paramsFinal <- combine_params_final(res)
#summary(paramsFinal)

# Boxplots grouped by nSims: one panel per parameter, one box per (stage, nSims).
plot_params_final(paramsFinal, par_thomas, type = type,
                  thinningType = thinningType, studyNr = studyNr, startParams = startParams,
                  title = "Parameter estimates by prefit stage and ensemble size")

# The same on the log scale the optimizer works on, where equal relative spread reads as
# equal visual spread -- easier to compare how fast each parameter tightens.
plot_params_final(paramsFinal, par_thomas, type = type,
                  thinningType = thinningType, studyNr = studyNr,
                  startParams = startParams, log_scale = TRUE, title = "As above, log scale")

# One parameter at full size.
plot_params_final(paramsFinal, par_thomas, type = type,
                  thinningType = thinningType, studyNr = studyNr,
                  startParams = startParams, parameters = "mu")

# How the spread falls with nSims at the end of the fit: ----
finalOnly <- paramsFinal[paramsFinal$prefitNumber == "final", ]
sapply(c("kappa", "omega", "mu"), function(pn)
  tapply(finalOnly[[pn]], finalOnly$nSims, sd))
