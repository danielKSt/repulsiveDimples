library(spatstat)
library(parallel)
library(mcprogress)
library(repulsiveDimples)

dataFolder <- "..."

# Functions to use in simulation study: ----
importance_sampling_thomas <- function(par_0, par_goal, repRange, nSims, xlims_sim, ylims_sim,
                                       nSamples, r_vec, printProgress = FALSE){
  # Unpack the parameter values
  kappa <- par_goal$kappa
  omega <- par_goal$omega
  mu <- par_goal$mu

  kappa_0 <- par_0$kappa
  omega_0 <- par_0$omega
  mu_0 <- par_0$mu

  if(printProgress){
    print("Simulating pattern for importance sampling:")
    patternSim <- pmclapply(X = rep(kappa_0, nSims), FUN = rThomas_matern_thinned,
                            scale = omega_0,
                            mu = mu_0,
                            repulsionRange = repRange,
                            xlims = xlims_sim,
                            ylims = ylims_sim,
                            saveparents = TRUE)
    print("Estimating baselines: ")
    K_lambda_baseline <- mcprogress::pmclapply(patternSim, estimate_K_lambda_baseline, r_vec = r_vec)
  } else {
    patternSim <- parallel::mclapply(X = rep(kappa_0, nSims), FUN = rThomas_matern_thinned,
                                     scale = omega_0,
                                     mu = mu_0,
                                     repulsionRange = repRange,
                                     win = simWin,
                                     saveparents = TRUE)
    K_lambda_baseline <- parallel::mclapply(patternSim, estimate_K_lambda_baseline, r_vec = r_vec)
  }

  rho_baseline <- sapply(patternSim, estimate_rho_baseline)
  if(printProgress){print("Calculating IS weights: ")}
  w_is <- importance_sampling_weigths(kappa = kappa, mu = mu, omega = omega,
                                      kappa_0 = kappa_0, mu_0 = mu_0, omega_0 = omega_0,
                                      patternSim = patternSim)$w_is

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
              K_is = K_is, K_is_norm = K_is_norm, w_is = w_is))
}

estimate_true_thomas <- function(par_goal, nTrue, xlims_true, ylims_true, repRange, r_vec, printProgress = FALSE){
  kappa <- par_goal$kappa
  omega <- par_goal$omega
  mu <- par_goal$mu

  startTime <- Sys.time()
  if(printProgress){
    print("Simulating with true parameters:")
    patternTrue <- pmclapply(X = rep(kappa, nTrue), FUN = rThomas_matern_thinned,
                             scale = omega,
                             mu = mu,
                             repulsionRange = repRange,
                             xlims = xlims_true,
                             ylims = ylims_true,
                             saveparents = TRUE)

    print("Estimating baselines: ")
    K_lambda_baseline <- mcprogress::pmclapply(patternTrue, estimate_K_lambda_baseline, r_vec = r_vec)
  } else {
    patternTrue <- parallel::mclapply(X = rep(kappa, nTrue), FUN = rThomas_matern_thinned,
                                      scale = omega,
                                      mu = mu,
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

  endTime <- Sys.time()
  return(list(rho_true = rho_true, K_true = K_true))
}

simulation_study_thomas_single_pars <- function(par_0, par_goal, repRange,
                                                nSims, xlims_sim, ylims_sim, nSamples, r_vec,
                                                res_true = NULL, nTrue = NULL,
                                                xlims_true = NULL, ylims_true = NULL,
                                                printProgress = FALSE){
  # Simulate with the true values to estimate
  if(is.null(res_true)){
    res_true <- estimate_true_thomas(par_goal = par_goal, nTrue = nTrue,
                                     xlims_true = xlims_true, ylims_true = ylims_true,
                                     repRange = repRange, printProgress = printProgress,
                                     r_vec = r_vec)
  }

  res_is <- importance_sampling_thomas(par_0 = par_0, par_goal = par_goal,
                                       repRange = repRange, nSims = nSims,
                                       xlims_sim = xlims_sim, ylims_sim = ylims_sim,
                                       nSamples = nSamples, r_vec = r_vec,
                                       printProgress = printProgress)


  return(list(res_is = res_is, res_true = res_true, nSamples = nSamples,
              par_0 = par_0, par_goal = par_goal, repRange = repRange,
              nSims = nSims, xlims_sim = xlims_sim, ylims_sim = ylims_sim,
              nTrue = nTrue, xlims_true = xlims_true, ylims_true = ylims_true))
}

combined_sim_study_thomas <- function(par_goal, par_delta, repRange,
                                      nSims, xlims_sim, ylims_sim, nSamples, r_vec,
                                      res_true = NULL, nTrue = NULL,
                                      xlims_true = NULL, ylims_true = NULL,
                                      printProgress = FALSE){
  if(is.null(res_true)){
    res <- vector(mode = "list", length = nrow(par_goal))
    for (i in 1:nrow(par_goal)) {
      print(paste("par_goal number:", i, " out of:", nrow(par_goal)))
      par_i <- par_goal[i, ]
      res_true <- estimate_true_thomas(par_goal = par_i, nTrue = nTrue, repRange = repRange,
                                       xlims_true = xlims_true, ylims_true = ylims_true,
                                       r_vec = r_vec, printProgress = printProgress)
      res_complete <- vector(mode = "list", length = nrow(par_delta))
      for(j in 1:nrow(par_delta)){
        print(paste("par_delta number:", j, " out of:", nrow(par_delta)))
        par_0 <- par_i + par_delta[j, ]

        res_complete[[j]] <- simulation_study_thomas_single_pars(par_0 = par_0, par_goal = par_i,
                                                                 repRange = repRange,
                                                                 nSims = nSims,
                                                                 xlims_sim = xlims_sim,
                                                                 ylims_sim = ylims_sim,
                                                                 nSamples = nSamples,
                                                                 r_vec = r_vec,
                                                                 res_true = res_true,
                                                                 printProgress = printProgress)
      }
      res[[i]] <- res_complete
    }
    return(res)
  }
}

# Run simulation study: ----

# Parameterar:

# Parameters and settings for simulation study 1:
repRange <- 0.4

par_goal <- data.frame(kappa = c(0.32),
                       omega = c(0.8),
                       mu = c(5.2))

par_delta <- data.frame(kappa = c(0.0, 0.01, 1.0),
                        omega = c(0.0, 0.0, 1.0),
                        mu = c(-0.8, 0.01, 1.0))

nSims <- 600000
xlims_sim <- c(0,6)
ylims_sim <- c(0,6)
nSamples <- floor(10^((seq(from = log10(10), to = log10(nSims), by = 0.333334))))
if(nSamples[length(nSamples)] < nSims){nSamples <- c(nSamples, nSims)}
r_vec <- c(seq(from = 0, to = 3, by = 0.1))

nTrue <- 200000
xlims_true <- c(0,6)
ylims_true <- c(0,6)

printProgress <- TRUE



# Parameters for currently running simulation study:

repRange <- 0.4

par_goal <- data.frame(kappa = c(0.32, 1.2),
                       omega = c(0.8, 1.2),
                       mu = c(5.2, 2.4))

par_delta <- data.frame(kappa = c(0.2, -0.2, 0.0, 0.0, 0.0, 0.0, 0.1),
                        omega = c(0.0, 0.0, 0.2, -0.2, 0.0, 0.0, 0.1),
                        mu = c(0.0, 0.0, 0.0, 0.0, 0.2, -0.2, 0.1))

nSims <- 60000
xlims_sim <- c(0,6)
ylims_sim <- c(0,6)
nSamples <- floor(10^((seq(from = log10(10), to = log10(nSims), by = 0.333334))))
if(nSamples[length(nSamples)] < nSims){nSamples <- c(nSamples, nSims)}
r_vec <- c(seq(from = 0, to = 3, by = 0.1))

nTrue <- 200000
xlims_true <- c(0,6)
ylims_true <- c(0,6)

printProgress <- TRUE

simStudyResults <- combined_sim_study_thomas(par_delta = par_delta,
                                             par_goal = par_goal, repRange = repRange,
                                             nSims = nSims, xlims_true = xlims_true, ylims_true = ylims_true,
                                             xlims_sim = xlims_sim, ylims_sim = ylims_sim, nSamples = nSamples,
                                             r_vec = r_vec, nTrue = nTrue,
                                             printProgress = TRUE)

save(simStudyResults, file = paste(dataFolder, "study2.RDa", sep = ""))
