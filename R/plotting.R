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


# save.pcf.plots <- function(re, we, bw_s, bw_d, l, spacing, rMax, dr_s, dr_d, wd, plot_folder){
#   if(is.null(wd)){
#     setwd("/Users/danielks/Library/CloudStorage/OneDrive-NTNU/PhD/Aalborg")
#   } else {
#     setwd(wd)
#   }
#
#   plot_folder <- "Figures/pcf/"
#
#   settings <- paste("RE", re, "_WE", we, sep = "")
#   load(file = paste("data/", settings, "/simulatedPoints.RDa", sep = ""))
#   load(file = paste("data/", settings, "/scales.RDa", sep = ""))
#
#   pcf_scars <- get_pcf(points_input = simulatedScars,
#                        l = l,
#                        spacing = spacing,
#                        bw = bw_s,
#                        rMax = rMax,
#                        dr = dr_s)
#
#   pcf_vorts <- get_pcf(points_input = simulatedVortices,
#                        l = l,
#                        spacing = spacing,
#                        bw = bw_d,
#                        rMax = rMax,
#                        dr = dr_d)
#
#   sPlot <- plot.pcf.ggplot(pcf_scars$g_r) + ggtitle(settings)
#   dPlot <- plot.pcf.ggplot(pcf_vorts$g_r) + ggtitle(settings)
#
#   ggsave(
#     filename = paste(plot_folder, settings, "vorts.pdf", sep = ""),
#     plot     = dPlot,
#     width    = 10, height = 8,
#     units    = "cm",
#     device   = cairo_pdf,
#     bg       = "transparent"
#   )
#   ggsave(
#     filename = paste(plot_folder, settings, "scars.pdf", sep = ""),
#     plot     = sPlot,
#     width    = 10, height = 8,
#     units    = "cm",
#     device   = cairo_pdf,
#     bg       = "transparent"
#   )
# }

# Kode for å lagre alle interressante pcf-plots

# save.pcf.plots(re = 2500, we = "inf",
#                bw_s = 0.08, bw_d = 0.04,
#                l = 256/lengthscale, spacing = 5, rMax = 2.5,
#                dr_s = 0.01, dr_d = 0.01)
# save.pcf.plots(re = 1000, we = "inf",
#                bw_s = 0.08, bw_d = 0.08,
#                l = 128/lengthscale, spacing = 5, rMax = 2.5,
#                dr_s = 0.01, dr_d = 0.01)
#
# save.pcf.plots(re = 2500, we = 10,
#                bw_s = 0.08, bw_d = 0.04,
#                l = 256/lengthscale, spacing = 5, rMax = 2.5,
#                dr_s = 0.01, dr_d = 0.01)
# save.pcf.plots(re = 1000, we = 10,
#                bw_s = 0.08, bw_d = 0.04,
#                l = 128/lengthscale, spacing = 5, rMax = 2.5,
#                dr_s = 0.01, dr_d = 0.01)
