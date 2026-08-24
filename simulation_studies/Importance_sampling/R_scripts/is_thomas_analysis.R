
library(ggplot2)
library(repulsiveDimples)
setwd("/Users/danielks/Library/CloudStorage/OneDrive-NTNU/PhD/Aalborg/repulsiveDimples/simulation_studies/")

load("study3a.RDa")

# Helper functions: ----
calc_K_diff <- function(K_true, case, q = 1/4){
  K_norm <- case$res_is$K_is_norm
  K_unnorm <- case$res_is$K_is
  dr <- K_true$r[2] - K_true$r[1]
  diff_norm <- c(1:length(K_norm))
  diff_unnorm <- c(1:length(K_unnorm))
  for(i in 1:length(K_norm)){
    diff_norm[i] <- sum((K_true$border^q - K_norm[[i]]$border^q)^2)*dr
    diff_unnorm[i] <- sum((K_true$border^q - K_unnorm[[i]]$border^q)^2)*dr
  }
  K_base <- case$res_is$K_base
  diff_base <- sum((K_true$border^q - K_base$border^q)^2)*dr
  return(list(normalized = diff_norm,
              unnormalized = diff_unnorm,
              base = diff_base))
}

# Basic study: ----

group_nr <- 2

group <- simStudyResults[[group_nr]]
res_true <- group$res_true

case_nr <- 6
case <- group$res_cases[[case_nr]]
print(exp(group$par_goal))
print(exp(case$par_0))

# Intensity
ymax <- max(c(res_true$rho_true, case$res_is$rho_is, case$res_is$rho_is_norm, case$res_is$rho_base))*1.1
ymin <- min(c(res_true$rho_true, case$res_is$rho_is, case$res_is$rho_is_norm, case$res_is$rho_base))*0.9

plot(x = log(group$nSamples), y = case$res_is$rho_is_norm, type = "l",
     ylim = c(ymin,ymax))
lines(x = log(group$nSamples), y = case$res_is$rho_is, lty = 2)
abline(h = res_true$rho_true, col = "red")
abline(h = case$res_is$rho_base, col = "purple")

# K-function
K_diff <- calc_K_diff(res_true$K_true, case)

print(K_diff$base)
plot(x = log(group$nSamples), y = K_diff$normalized, type = "l")
lines(x = log(group$nSamples), y = K_diff$unnormalized, lty = 2)
abline(h = K_diff$base, col = "red")

plot(x = res_true$K_true$r, y = res_true$K_true$border^0.25, type = "l")
lines(x = res_true$K_true$r, y = case$res_is$K_base$border^0.25, col = "red")
lines(x = res_true$K_true$r, y = case$res_is$K_is[[13]]$border^0.25, col = "purple", lty = 2)
lines(x = res_true$K_true$r, y = case$res_is$K_is_norm[[13]]$border^0.25, col = "orange", lty = 2)

# Look at effective sample size
ess <- 1/sum((case$res_is$w_is/sum(case$res_is$w_is))^2)

ess_vec <- sapply(X = group$res_cases, FUN = function(case) 1/sum((case$res_is$w_is/sum(case$res_is$w_is))^2))
ess_vec[2:length(ess_vec)]
which.max(ess_vec[2:length(ess_vec)])

exp(group$res_cases[[which.max(ess_vec[2:length(ess_vec)])+1]]$par_0)
exp(group$par_goal)

# Window study: ----
rm(list = ls())
load("studyWindow.RDa")

window_case <- simStudyResults[[1]]
calc_ess <- function(window_case){
  window_case <- window_case[[1]]
  ess_vec <- sapply(X = window_case$res_cases, FUN = function(case) 1/sum((case$res_is$w_is/sum(case$res_is$w_is))^2))
  return(ess_vec)
}

ess <- lapply(X = simStudyResults, calc_ess)

ess[[1]]
ess[[2]]
ess[[3]]
ess[[4]]
ess[[5]]
ess[[6]]

# Kappa study: ----
rm(list = ls())
load("studyKappa.RDa")

calc_ess <- function(window_case){
  #window_case <- window_case[[1]]
  ess_vec <- sapply(X = window_case$res_cases, FUN = function(case) 1/sum((case$res_is$w_is/sum(case$res_is$w_is))^2))
  return(ess_vec)
}
ess <- lapply(X = simStudyResults, calc_ess)

ess_max <- sapply(ess, max)
ess_max
