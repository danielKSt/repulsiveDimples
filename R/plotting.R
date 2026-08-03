#' Plot pcf using ggplot2
#'
#' @description
#' Plot a pcf with ggplot, ignoring the first r-value in case this is NA
#'
#' @param pcf_input Pcf function to be plotted
#'
#' @importFrom rlang .data
#'
#' @export
mypcfplot.ggplot <- function(pcf_input){
  pcf_input <- pcf_input[2:nrow(pcf_input), ]
  df <- data.frame(r = pcf_input$r, y = pcf_input$theo, type = "theo")
  df <- rbind(df, data.frame(r = pcf_input$r, y = pcf_input$iso, type = "iso"))
  df <- rbind(df, data.frame(r = pcf_input$r, y = pcf_input$trans, type = "trans"))
  res <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = .data$r, y = .data$y, colour = .data$type, linetype = .data$type))+
    ggplot2::geom_line()+
    ggplot2::scale_color_manual(values = c(iso = "steelblue", trans = "firebrick", theo = "darkgreen"))+
    ggplot2::scale_linetype_manual(values = c(iso = "solid", trans = "solid", theo = "dashed"))
  return(res)
}


#' Compare the results from optimization with the target K-function and intensity
#'
#' @description
#' Compare the results from optimization with the target K-function and intensity
#' @param rho_hat Intensity estimated from the data
#' @param K_hat K-function estimated from the data
#' @param param_estim The parameters estimated from the minimum contrast estimation
#' @param repRange_estim Estimated repulsion range
#' @param q Power to apply to K-function
#' @param nSims numer of simulations
#' @param xrange xrange for simulation window
#' @param yrange yrange for simulation window
#'
#' @export
result_verification <- function(rho_hat, K_hat, param_estim, repRange_estim, q = 1/4,
                                nSims = 1000, xrange = c(0,10), yrange = c(0,10)){
  verifPattern <- lapply(rep(param_estim[1], nSims), rThomas_matern_thinned,
                         scale = param_estim[2], mu = param_estim[3],
                         repulsionRange = repRange_estim, xlims = xrange, ylims = yrange)

  rho_verif_baseline <- sapply(verifPattern, estimate_rho_baseline)
  K_lambda_verif_baseline <- lapply(verifPattern, estimate_K_lambda_baseline, r_vec = K_hat$r)

  w_is <- rep(1, length(verifPattern))

  rho_verif <- rho_importance_sampling(w_is = rep(1, length(verifPattern)), rho_baseline = rho_verif_baseline)
  K_verif <- K_importance_sampling(w_is = rep(1, length(verifPattern)),
                                   K_lambda_baseline = K_lambda_verif_baseline,
                                   rho_baseline = rho_verif_baseline)

  print(paste("Target intensity: ", rho_hat))
  print(paste("Parameter intensity: ", rho_verif))

  graphics::plot(x = K_hat$r, y = K_hat$border^q, type = "l", ylim = c(0,max(c(K_hat$border^q, K_verif$border^q))),
       ylab = "K(r)", xlab = "r")
  graphics::lines(x = K_verif$r, y = K_verif$border^q, col = "orange", lty = 2)
  graphics::legend("topleft", c("Target K-function", "Parameter K-function"), col = c("black", "orange"), lty = c(1, 2))
}
