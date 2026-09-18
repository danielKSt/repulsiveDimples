
library(ggplot2)
library(viridis)
library(dplyr)
library(patchwork)
library(repulsiveDimples)
# Run from the project root. The paths below are relative to this directory, so a second
# sourcing in the same session fails here, the working directory already being it.
setwd("simulation_studies/Importance_sampling/")

source("R_scripts/is_analysis_helpers.R")
load("Results/study4.RDa")

# Simple convergence plots: ----

group <- simStudyResults[[1]]
kappas <- sapply(group$res_cases, function(case) case$par_0$kappa)
kappasLogDiff <- group$par_goal$kappa - kappas
cases <- intersect(which(abs(kappasLogDiff) > 0.001), which(abs(kappasLogDiff) < 0.3))


res <- df_convergence_results_complete(group = group, basecase_ind = 1, nSamplesMin = 200)

res_kappa <- res |>
  filter(abs(kappasLogDiff) > 0.001) |>
  filter(abs(kappasLogDiff) < 0.3)


rho_norm <- ggplot(data = res_kappa, aes(x = log(nSamples), y = rho_is_norm, group = case_ind, colour = kappasLogDiff)) +
  geom_line() +
  geom_hline(yintercept = group$res_true$rho_true, linetype = "dashed", colour = "red") +
  scale_color_viridis_c(name = bquote("log difference in " * kappa), option = "H") +
  ylab(bquote(hat(rho)[Norm])) +
  ylim(c(0.9, 1.15))

rho_unnorm <- ggplot(data = res_kappa, aes(x = log(nSamples), y = rho_is, group = case_ind, colour = kappasLogDiff)) +
  geom_line() +
  geom_hline(yintercept = group$res_true$rho_true, linetype = "dashed", colour = "red") +
  scale_color_viridis_c(name = bquote("log difference in " * kappa), option = "H") +
  ylab(bquote(hat(rho)[Unnorm])) +
  ggtitle(bquote("Estimates for " * rho * ", with " *
                   kappa * " = " * .(round(exp(group$par_goal$kappa), digits = 6)) * ", "  *
                   sigma^2 * " = " * .(round(exp(group$par_goal$omega), digits = 6)) * ", " *
                   mu * " = " * .(round(exp(group$par_goal$mu), digits = 6)))) +
  ylim(c(0.9, 1.15))

rho_plot <- rho_unnorm + rho_norm + plot_layout(guides = "collect")

K_unnorm <- ggplot(data = res_kappa, aes(x = log(nSamples), y = K_diff, group = case_ind, colour = kappasLogDiff)) +
  geom_line() +
  scale_color_viridis_c(name = bquote("log difference in " * kappa), option = "H") +
  ylab(bquote(hat(K)[Unnorm])) +
  ylim(c(0,0.007)) +
  ggtitle(bquote("Estimates for " *
                   kappa * " = " * .(round(exp(group$par_goal$kappa), digits = 6)) * ", "  *
                   sigma^2 * " = " * .(round(exp(group$par_goal$omega), digits = 6)) * ", " *
                   mu * " = " * .(round(exp(group$par_goal$mu), digits = 6))))

K_norm <- ggplot(data = res_kappa, aes(x = log(nSamples), y = K_diff_norm, group = case_ind, colour = kappasLogDiff)) +
  geom_line() +
  scale_color_viridis_c(name = bquote("log difference in " * kappa), option = "H") +
  ylab(bquote(hat(K)[Norm])) +
  ylim(c(0,0.007))

K_plot <- K_unnorm + K_norm + plot_layout(guides = "collect")

plotFolder <- "/Users/danielks/Library/CloudStorage/OneDrive-NTNU/PhD/Aalborg/presentasjonar/figs/is_sim_studies/"
ggsave(
  filename = paste(plotFolder, "study4_group1_rho.pdf", sep = ""),
  plot     = rho_plot,
  width    = 32, height = 16,
  units    = "cm",
  device   = cairo_pdf,
  bg       = "transparent"
)

ggsave(
  filename = paste(plotFolder, "study4_group1_K.pdf", sep = ""),
  plot     = K_plot,
  width    = 32, height = 16,
  units    = "cm",
  device   = cairo_pdf,
  bg       = "transparent"
)


res_sigma <- res |>
  filter(abs(sigmasqsLogDiff) > 0.001) |>
  filter(abs(sigmasqsLogDiff) < 0.3)

ggplot(data = res_sigma, aes(x = log(nSamples), y = rho_is_norm, group = case_ind, colour = sigmasqsLogDiff)) +
  geom_line() +
  geom_hline(yintercept = group$res_true$rho_true, linetype = "dashed", colour = "red") +
  scale_color_viridis_c(name = bquote("log difference in " * kappa), option = "H") +
  ylab(bquote(hat(rho)[Norm])) +
  ylim(c(0.9, 1.15))

res_mu <- res |>
  filter(abs(musLogDiff) > 0.001) |>
  filter(abs(musLogDiff) < 0.3)

ggplot(data = res_mu, aes(x = log(nSamples), y = rho_is_norm, group = case_ind, colour = musLogDiff)) +
  geom_line() +
  geom_hline(yintercept = group$res_true$rho_true, linetype = "dashed", colour = "red") +
  scale_color_viridis_c(name = bquote("log difference in " * kappa), option = "H") +
  ylab(bquote(hat(rho)[Norm])) +
  ylim(c(0.9, 1.15))

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
