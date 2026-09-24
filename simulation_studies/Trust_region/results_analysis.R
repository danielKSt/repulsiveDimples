library(spatstat)
library(ggplot2)
library(repulsiveDimples)

source(file = "simulation_studies/Trust_region/Helper_scripts/results_analysis_helpers.R")

studyNr <- 1

# Fetch start values for the parameters: ----
load(file = paste0("simulation_studies/Trust_region/Data/study", studyNr, ".RDa"))

start_fitted <- thomas.estK(X = K_hat_unions, rmin = 2*par_thomas$rRange, rmax = 3)
lambda <- -log(1-pi*rho_hat*par_thomas$rRange^2)/(pi*par_thomas$rRange^2)
start_fitted <- start_fitted$par

start_params <- log(c(start_fitted[1], start_fitted[2], lambda/start_fitted[1]))
names(start_params) <- c("kappa", "omega", "mu")
rm(start_fitted, lambda, K_hat_unions, nTrue, sidelength, nUnions)

# Look at the results for a single run: ----
load(file = paste0("simulation_studies/Trust_region/Results/Single_run/study", studyNr, "_res.RDa"))

a <- getParamsFinal(res = res, nSims = 40000)
a

ggplot(a, aes(prefitNumber, kappa)) + geom_boxplot() + geom_hline(yintercept = par_thomas$kappa, colour = "red")
ggplot(a, aes(prefitNumber, omega)) + geom_boxplot() + geom_hline(yintercept = par_thomas$omega, colour = "red")
ggplot(a, aes(prefitNumber, mu)) + geom_boxplot() + geom_hline(yintercept = par_thomas$mu, colour = "red")

rm(res, a)

# Look at the results with multiple runs: ----
load(file = paste0("simulation_studies/Trust_region/Results/study", studyNr, "_main.RDa"))

# study1_main.RDa predates run_level returning list(res, nSims), so it holds res1..res4
# and has to be assembled into the shape combineParamsFinal expects. Files written by
# study2_run.R onwards already store `res` in that shape, and can be passed straight in.
if(studyNr == 1){
  res <- list(list(res = res1, nSims = 500),
              list(res = res2, nSims = 1000),
              list(res = res3, nSims = 5000),
              list(res = res4, nSims = 10000))
  rm(res1, res2, res3, res4)
}
paramsFinal <- combineParamsFinal(res)
summary(paramsFinal)

# Boxplots grouped by nSims: one panel per parameter, one box per (stage, nSims).
plotParamsFinal(paramsFinal, par_thomas, start_params) +
  labs(title = "Parameter estimates by prefit stage and ensemble size")

# The same on the log scale the optimizer works on, where equal relative spread reads as
# equal visual spread -- easier to compare how fast each parameter tightens.
plotParamsFinal(paramsFinal, par_thomas, start_params, log_scale = TRUE) +
  labs(title = "As above, log scale")

# One parameter at full size.
plotParamsFinal(paramsFinal, par_thomas, start_params, parameters = "kappa")

# How the spread falls with nSims at the end of the fit: ----
finalOnly <- paramsFinal[paramsFinal$prefitNumber == "final", ]
sapply(c("kappa", "omega", "mu"), function(pn)
  tapply(finalOnly[[pn]], finalOnly$nSims, sd))
