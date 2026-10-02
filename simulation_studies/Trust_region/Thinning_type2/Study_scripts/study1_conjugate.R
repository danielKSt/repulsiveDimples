# study1, fitted with conjugate direction line searches (Powell 1964).
#
# It is paired with the quadratic fits in Results/study1_quad.RDa, which quad_run.R makes
# when type 2 is included. Those run the same data, the same starting parameters and the
# same seeds through the other trust step method, so the two sets of fits differ only in
# how each step searches. Results go to Results/study1_conjugate.RDa, with one file per
# ensemble size alongside it.

# Load data and libraries: ----
library(spatstat)
library(parallel)
library(mcprogress)
library(repulsiveDimples)

source("simulation_studies/Trust_region/Helper_scripts/study_run_helpers.R")

# Matern thinning type to fit. It also decides which Thinning_type<N> folder the data are
# loaded from and the results saved in.
thinningType <- 2
load(file = paste0(study_folder(thinningType), "Data/study1.RDa"))

# Run the optimizer: ----
start_params <- study_start_params(K_hat_unions = K_hat_unions, par_thomas = par_thomas,
                                   rho_hat = rho_hat)
rm(K_hat_unions, nTrue, sidelength, nUnions)

res <- run_study(study = "study1", variant = "conjugate", method = "conjugate",
                 start_params = start_params, par_thomas = par_thomas,
                 rho_hat = rho_hat, K_hat = K_hat, thinningType = thinningType,
                 max.iter.prefit = 10)
