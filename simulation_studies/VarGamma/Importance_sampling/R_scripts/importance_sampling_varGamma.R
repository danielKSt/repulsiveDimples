library(spatstat)
library(parallel)
library(mcprogress)
library(repulsiveDimples)

# mclapply() defaults to mc.cores = 2 regardless of how many cores the machine
# has, which silently caps every parallel step in these studies at 2x. Set it
# once here so all sourcing run-scripts inherit it. Set options(mc.cores = ...)
# before sourcing this file to override (e.g. to a scheduler's allocation
# rather than the whole machine).
if (is.null(getOption("mc.cores"))) {
  options(mc.cores = parallel::detectCores())
}

# Under RNGkind("L'Ecuyer-CMRG"), which the run scripts use so that set.seed makes them
# reproducible, mclapply hands its workers streams derived from the main process's seed,
# and nothing in the call advances that seed. Two simulation calls in a row would then
# reuse the same random numbers, making the ensemble for every case a copy of the start
# of the one behind the true values. Drawing a fresh seed in the main process before each
# simulation keeps the ensembles independent, and the run still reproducible.
reseed_main_process <- function(){
  set.seed(sample.int(.Machine$integer.max, 1))
}

# Functions to use in simulation study: ----
# The Variance gamma counterparts of importance_sampling_thomas.R (in
# simulation_studies/Thomas/Importance_sampling/R_scripts), with the same
# structure and the same returned lists, so the two studies can be analysed side by
# side. The parameters are (kappa, scale, mu) on the log scale, in data frames with
# those column names. The shape nu is fixed: it is passed separately, never perturbed,
# and stored with the results.
importance_sampling_varGamma <- function(par_0, par_goal, nu, repRange, nSims, xlims_sim, ylims_sim,
                                         nSamples, r_vec, printProgress = FALSE){
  # Unpack parameters and transform them out of log scale
  kappa <- exp(par_goal$kappa)
  scale <- exp(par_goal$scale)
  mu <- exp(par_goal$mu)

  kappa_0 <- exp(par_0$kappa)
  scale_0 <- exp(par_0$scale)
  mu_0 <- exp(par_0$mu)

  reseed_main_process()
  if(printProgress){
    print("Simulating pattern for importance sampling:")
    patternSim <- pmclapply(X = rep(kappa_0, nSims), FUN = rVarGamma_matern_thinned,
                            scale = scale_0,
                            mu = mu_0,
                            nu = nu,
                            repulsionRange = repRange,
                            xlims = xlims_sim,
                            ylims = ylims_sim,
                            saveparents = TRUE)
    print("Estimating baselines: ")
    K_lambda_baseline <- mcprogress::pmclapply(patternSim, estimate_K_lambda_baseline, r_vec = r_vec)
  } else {
    patternSim <- parallel::mclapply(X = rep(kappa_0, nSims), FUN = rVarGamma_matern_thinned,
                                     scale = scale_0,
                                     mu = mu_0,
                                     nu = nu,
                                     repulsionRange = repRange,
                                     xlims = xlims_sim,
                                     ylims = ylims_sim,
                                     saveparents = TRUE)
    K_lambda_baseline <- parallel::mclapply(patternSim, estimate_K_lambda_baseline, r_vec = r_vec)
  }

  rho_baseline <- sapply(patternSim, estimate_rho_baseline)
  if(printProgress){print("Calculating IS weights: ")}
  w_is <- importance_sampling_weights_varGamma(kappa = kappa, mu = mu, scale = scale, nu = nu,
                                               kappa_0 = kappa_0, mu_0 = mu_0, scale_0 = scale_0,
                                               patternSim = patternSim, parallel = TRUE)$w_is

  rho_base <- rho_importance_sampling(rho_baseline = rho_baseline, w_is = rep(1, nSims))
  K_base <- K_importance_sampling(w_is = rep(1, nSims), rho_baseline = rho_baseline,
                                  K_lambda_baseline = K_lambda_baseline, r_vec = r_vec)

  rho_is <- c(1:length(nSamples))
  rho_is_norm <- c(1:length(nSamples))
  K_is <- vector(mode = "list", length = length(nSamples))
  K_is_norm <- vector(mode = "list", length = length(nSamples))

  for(i in 1:length(nSamples)){
    # Estimate intensity:
    rho_is[i] <- rho_importance_sampling(rho_baseline = rho_baseline[1:nSamples[i]],
                                         w_is = w_is[1:nSamples[i]])

    rho_is_norm[i] <- rho_importance_sampling(rho_baseline = rho_baseline[1:nSamples[i]],
                                              w_is = w_is[1:nSamples[i]],
                                              normalized = TRUE)

    # Estimate K-function
    K_is[[i]] <- K_importance_sampling(w_is = w_is[1:nSamples[i]], rho_baseline = rho_baseline[1:nSamples[i]],
                                       K_lambda_baseline = K_lambda_baseline[1:nSamples[i]], r_vec = r_vec)

    K_is_norm[[i]] <- K_importance_sampling(w_is = w_is[1:nSamples[i]], rho_baseline = rho_baseline[1:nSamples[i]],
                                            K_lambda_baseline = K_lambda_baseline[1:nSamples[i]], r_vec = r_vec, normalized = TRUE)
  }
  return(list(rho_is = rho_is, rho_is_norm = rho_is_norm,
              K_is = K_is, K_is_norm = K_is_norm, w_is = w_is,
              rho_base = rho_base, K_base = K_base))
}

estimate_true_varGamma <- function(par_goal, nu, nTrue, xlims_true, ylims_true, repRange, r_vec, printProgress = FALSE){
  # Unpack parameters and transform them out of log scale
  kappa <- exp(par_goal$kappa)
  scale <- exp(par_goal$scale)
  mu <- exp(par_goal$mu)

  reseed_main_process()
  if(printProgress){
    print("Simulating with true parameters:")
    patternTrue <- pmclapply(X = rep(kappa, nTrue), FUN = rVarGamma_matern_thinned,
                             scale = scale,
                             mu = mu,
                             nu = nu,
                             repulsionRange = repRange,
                             xlims = xlims_true,
                             ylims = ylims_true,
                             saveparents = TRUE)

    print("Estimating baselines: ")
    K_lambda_baseline <- mcprogress::pmclapply(patternTrue, estimate_K_lambda_baseline, r_vec = r_vec)
  } else {
    patternTrue <- parallel::mclapply(X = rep(kappa, nTrue), FUN = rVarGamma_matern_thinned,
                                      scale = scale,
                                      mu = mu,
                                      nu = nu,
                                      repulsionRange = repRange,
                                      xlims = xlims_true,
                                      ylims = ylims_true,
                                      saveparents = TRUE)
    K_lambda_baseline <- parallel::mclapply(patternTrue, estimate_K_lambda_baseline, r_vec = r_vec)
  }

  rho_baseline <- sapply(patternTrue, estimate_rho_baseline)

  if(printProgress){print("Estimating rho:")}
  rho_true <- mean(rho_baseline)
  if(printProgress){print("Estimating K:")}
  K_true <- K_importance_sampling(w_is = rep(1, nTrue), K_lambda_baseline = K_lambda_baseline,
                                  rho_baseline = rho_baseline, normalized = TRUE)

  return(list(rho_true = rho_true, K_true = K_true, nTrue = nTrue,
              xlims_true = xlims_true, ylims_true = ylims_true))
}


combined_sim_study_varGamma <- function(par_goal, par_delta, nu, repRange,
                                        nSims, xlims_sim, ylims_sim, nSamples, r_vec,
                                        nTrue, xlims_true, ylims_true,
                                        printProgress = FALSE){
  res <- vector(mode = "list", length = nrow(par_goal))
  for (i in 1:nrow(par_goal)) {
    print(paste("par_goal number:", i, " out of:", nrow(par_goal)))
    par_i <- par_goal[i, ]
    res_true <- estimate_true_varGamma(par_goal = par_i, nu = nu, nTrue = nTrue, repRange = repRange,
                                       xlims_true = xlims_true, ylims_true = ylims_true,
                                       r_vec = r_vec, printProgress = printProgress)
    res_complete <- vector(mode = "list", length = nrow(par_delta))
    for(j in 1:nrow(par_delta)){
      print(paste("par_delta number:", j, " out of:", nrow(par_delta)))
      par_0 <- par_i + par_delta[j, ]

      res_is <- importance_sampling_varGamma(par_0 = par_0, par_goal = par_i, nu = nu,
                                             repRange = repRange, nSims = nSims,
                                             xlims_sim = xlims_sim, ylims_sim = ylims_sim,
                                             nSamples = nSamples, r_vec = r_vec,
                                             printProgress = printProgress)

      res_complete[[j]] <- list(res_is = res_is, par_0 = par_0, nu = nu, repRange = repRange,
                                nSims = nSims, xlims_sim = xlims_sim, ylims_sim = ylims_sim)
    }
    res[[i]] <- list(res_true = res_true,
                     nSamples = nSamples,
                     par_goal = par_i,
                     nu = nu,
                     res_cases = res_complete)
  }
  return(res)
}
