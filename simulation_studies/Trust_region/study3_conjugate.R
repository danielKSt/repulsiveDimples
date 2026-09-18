# study3, fitted with conjugate direction line searches (Powell 1964).
#
# The paired script is study3_quad.R, which runs the same data, the same starting
# parameters and the same seeds through the other trust step method, so the two sets of
# fits differ only in how each step searches. Results go to Results/study3_conjugate.RDa,
# with one file per ensemble size alongside it.

# Load data and libraries: ----
library(spatstat)
library(parallel)
library(mcprogress)
library(repulsiveDimples)

source("simulation_studies/Trust_region/study_run_helpers.R")
dir.create(resFolder, showWarnings = FALSE, recursive = TRUE)

load(file = "simulation_studies/Trust_region/Data/study3.RDa")

# Run the optimizer: ----
start_params <- study_start_params(K_hat_unions = K_hat_unions, par_thomas = par_thomas,
                                   rho_hat = rho_hat)
rm(K_hat_unions, nTrue, sidelength, nUnions)

res <- run_study(study = "study3", variant = "conjugate", method = "conjugate",
                 start_params = start_params, par_thomas = par_thomas,
                 rho_hat = rho_hat, K_hat = K_hat,
                 max.iter.prefit = 15)
