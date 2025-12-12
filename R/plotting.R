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

