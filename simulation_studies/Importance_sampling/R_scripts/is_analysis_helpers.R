

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


df_convergence_results <- function(group, cases, basecase_ind, paramLogDiff, nSamplesMin){
  K_diff <- calc_K_diff(K_true = group$res_true$K_true, case = group$res_cases[[basecase_ind]])
  inds <- which(group$nSamples > nSamplesMin)
  res <- data.frame(nSamples = group$nSamples[inds],
                    rho_is_norm = group$res_cases[[basecase_ind]]$res_is$rho_is_norm[inds],
                    rho_is = group$res_cases[[basecase_ind]]$res_is$rho_is[inds],
                    K_diff_norm = K_diff$normalized[inds],
                    K_diff = K_diff$unnormalized[inds],
                    case_ind = basecase_ind,
                    paramLogDiff = paramLogDiff[basecase_ind])
  for (case in cases) {
    K_diff <- calc_K_diff(K_true = group$res_true$K_true, case = group$res_cases[[case]])
    res <- rbind(res,
                 data.frame(nSamples = group$nSamples[inds],
                            rho_is_norm = group$res_cases[[case]]$res_is$rho_is_norm[inds],
                            rho_is = group$res_cases[[case]]$res_is$rho_is[inds],
                            K_diff_norm = K_diff$normalized[inds],
                            K_diff = K_diff$unnormalized[inds],
                            case_ind = case,
                            paramLogDiff = paramLogDiff[case]))
  }
  return(res)
}

df_convergence_results_complete <- function(group, cases = NULL, basecase_ind, nSamplesMin){
  # Preparations:
  if(is.null(cases)){
    cases <- which(1:length(group$res_cases) != basecase_ind)
  }

  kappas <- sapply(group$res_cases, function(case) case$par_0$kappa)
  kappasLogDiff <- group$par_goal$kappa - kappas

  sigmasqs <- sapply(group$res_cases, function(case) case$par_0$omega)
  sigmasqsLogDiff <- group$par_goal$omega - sigmasqs

  mus <- sapply(group$res_cases, function(case) case$par_0$mu)
  musLogDiff <- group$par_goal$mu - mus

  inds <- which(group$nSamples > nSamplesMin)

  # Basecase:
  K_diff <- calc_K_diff(K_true = group$res_true$K_true, case = group$res_cases[[basecase_ind]])
  res <- data.frame(nSamples = group$nSamples[inds],
                    rho_is_norm = group$res_cases[[basecase_ind]]$res_is$rho_is_norm[inds],
                    rho_is = group$res_cases[[basecase_ind]]$res_is$rho_is[inds],
                    K_diff_norm = K_diff$normalized[inds],
                    K_diff = K_diff$unnormalized[inds],
                    case_ind = basecase_ind,
                    kappasLogDiff = kappasLogDiff[basecase_ind],
                    sigmasqsLogDiff = sigmasqsLogDiff[basecase_ind],
                    musLogDiff = musLogDiff[basecase_ind])
  for (case in cases) {
    K_diff <- calc_K_diff(K_true = group$res_true$K_true, case = group$res_cases[[case]])
    res <- rbind(res,
                 data.frame(nSamples = group$nSamples[inds],
                            rho_is_norm = group$res_cases[[case]]$res_is$rho_is_norm[inds],
                            rho_is = group$res_cases[[case]]$res_is$rho_is[inds],
                            K_diff_norm = K_diff$normalized[inds],
                            K_diff = K_diff$unnormalized[inds],
                            case_ind = case,
                            kappasLogDiff = kappasLogDiff[case],
                            sigmasqsLogDiff = sigmasqsLogDiff[case],
                            musLogDiff = musLogDiff[case]))
  }
  return(res)
}
