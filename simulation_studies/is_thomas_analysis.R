
library(ggplot2)
setwd("/Users/danielks/Library/CloudStorage/OneDrive-NTNU/PhD/Aalborg/repulsiveDimples/simulation_studies/")

load("study2.RDa")

# Helper functions: ----
calc_K_diff <- function(K_true, K_norm, K_unnorm, q = 1/4){
  dr <- K_true$r[2] - K_true$r[1]
  diff_norm <- c(1:length(K_norm))
  diff_unnorm <- c(1:length(K_unnorm))
  for(i in 1:length(K_norm)){
    diff_norm[i] <- sum((K_true$border^q - K_norm[[i]]$border^q)^2)*dr
    diff_unnorm[i] <- sum((K_true$border^q - K_unnorm[[i]]$border^q)^2)*dr
  }
  return(list(normalized = diff_norm,
              unnormalized = diff_unnorm))
}


# Overview: ----
for(group_nr in 1:length(simStudyResults)){
  group <- simStudyResults[[group_nr]]
  print("#######################################################################")
  print(paste("----------- Analysis of results from group", group_nr, ": -----------"))
  print("#######################################################################")
  res_true <- group[[1]]$res_true
  par_goal <- group[[1]]$par_goal

  print("Some information from this group:")
  print(paste("kappa_goal:", par_goal$kappa))
  print(paste("omega_goal:", par_goal$omega))
  print(paste("mu_goal:", par_goal$mu))
  print(paste("rho_true:", res_true$rho_true))


  for(case_nr in 1:length(group)){
    case <- group[[case_nr]]
    par_0 <- case$par_0

    print("")
    print(paste("----------- Analysis of results from case", case_nr, ": -----------"))

    print("Some information from this case:")
    print(paste("kappa_0:", par_goal$kappa))
    print(paste("omega_0:", par_goal$omega))
    print(paste("mu_0:", par_goal$mu))

    # print(paste("rho:", res_true$rho_true)) not available in study2
    rho_diff <- res_true$rho_true - case$res_is$rho_is
    rho_diff_norm <- res_true$rho_true - case$res_is$rho_is_norm
    print("Differences in rho with unnormalized IS weights:")
    print(rho_diff)
    print("Differences in rho with normalized IS weights:")
    print(rho_diff_norm)

    K_diff <- calc_K_diff(K_true = res_true$K_true, K_norm = case$res_is$K_is_norm, K_unnorm = case$res_is$K_is)

    print("Contrast of K-function with unnormalized IS weights:")
    print(K_diff$unnormalized)
    print("Contrast of K-function with normalized IS weights:")
    print(K_diff$normalized)
  }
  print("")
  print("")
}

# Visual inspections

group_nr <- 1
case_nr <- 2

group <- simStudyResults[[group_nr]]
case <- simStudyResults[[group_nr]][[1]]

rho_df <- data.frame(x = log(case$nSamples), rho_unnorm = case$res_is$rho_is,
                     rho_norm = case$res_is$rho_is_norm, case = 1)
K_df <- data.frame(x = log(nSamples), )
for(i in 2:length(group)){
  rho_df <- rbind(rho_df,
                  data.frame(x = log(group[[i]]$nSamples), rho_unnorm = group[[i]]$res_is$rho_is,
                             rho_norm = group[[i]]$res_is$rho_is_norm, case = i))
}

ggplot(data = rho_df, mapping = aes(x = x, y = rho_unnorm, colour = case)) +
  geom_line() +
  #geom_line(mapping = aes(x = x, y = rho_norm, colour = case), linetype = "dashed") +
  geom_hline(yintercept = case$res_true$rho_true)

plot(x = log(case$nSamples), y = case$res_is$rho_is_norm, type = "l")
lines(x = log(case$nSamples), y = case$res_is$rho_is, lty = 2)
abline(h = case$res_true$rho_true, col = "red")
plot(x = log(case$nSamples), y = case$res_is$rho_is_norm, type = "l")

