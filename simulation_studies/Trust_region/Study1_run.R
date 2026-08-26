# Load data and libraries: ----
library(spatstat)
library(ggplot2)
library(parallel)
library(mcprogress)
devtools::load_all()

load(file = "simulation_studies/Trust_region/Data/study1.RDa")

# Run the optimizer: ----

# Pre-fit as if unthinned:
start_fitted <- thomas.estK(X = K_hat_unions, rmin = 2*min_dist, rmax = 3)
lambda <- -log(1-pi*rho_hat*min_dist^2)/(pi*min_dist^2)
start_fitted <- start_fitted$par

start_params <- log(c(start_fitted[1], start_fitted[2], lambda/start_fitted[1]))
rm(start_fitted, lambda, K_hat_unions, nTrue, sidelength)

#params = c(log(start_params[1]), log(par_thomas$omega), log(par_thomas$mu))

res <- min_contrast_trust_region(params = start_params,
                                 parFreeIndex = c(1, 2, 3), repRange = min_dist, rho_hat = rho_hat,
                                 K_hat = K_hat, xlims = c(0,3), ylims = c(0,3), nSims = 40000, deltaInit = 0.2,
                                 eta = 0.05, deltaMax = 1.0, wq = c(1000, 1/4), normalized = TRUE, deltaMin = 0.00001,
                                 tol = 10^-8, max.iter = 50, max.iter.prefit = 50, printProgress = TRUE)


par_res <- exp(c(res$params[nrow(res$params), 1], log(par_thomas$omega), log(par_thomas$mu)))
par_thomas
start_params
par_res


# Assign parameters manually for hand-testing:
params <- start_params
parFreeIndex <- c(1, 2, 3)
repRange <- min_dist
xlims <- c(0,3)
ylims <- c(0,3)
nSims <- 40000
deltaInit <- 0.2
deltaMax <- 1.0
deltaMin <- 10^-6
wq <- c(1000, 1/4)
normalized <- TRUE
tol <- 10^-8
max.iter <- 50
max.iter.prefit <- 50
printProgress <- TRUE
eta <- 0.05
