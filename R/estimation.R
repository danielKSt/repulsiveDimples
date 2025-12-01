
#' Function for estimating the intensity given a list of point patterns
#' The point patterns are assumed to come from a time series of point patterns
#'
#' @param points_input List of point patterns used for estimation
#' @param l Side length of observation window
#' @param timescale How many indices apart do snapshots need to be in order to be independent?
#'
#' @export
intensity_est <- function(points_input, l, timescale){
  snapshots <- seq(from = 1, to = length(points_input), by = ceiling(timescale))
  points.ppp <- vector(mode = "list", length = length(snapshots))
  nPoints <- 0

  for (t_ind in 1:length(snapshots)) {
    t <- snapshots[t_ind]
    if(!is.numeric(points_input[[t]])){
      nPoints <- nPoints + length(points_input[[t]]$x)
    }
  }
  l_hat <- nPoints/(length(snapshots)*l^2)

  return(l_hat)
}

#' Function for estimating the K-function given a list of point patterns
#' The point patterns are assumed to come from a time series of point patterns
#'
#' @param points_input List of point patterns used for estimation
#' @param l Side length of observation window
#' @param spacing How far apart are observations window in order to be independent?
#' @param rMax Max radius for which estimation is performed
#' @param dr How fine grid for estimation?
#' @param timescale How many indices apart do snapshots need to be in order to be independent?
#'
#' @export
K_est.unions <- function(points_input, l, spacing, rMax, dr, timescale){
  snapshots <- seq(from = 1, to = length(points_input), by = ceiling(timescale))
  points.ppp <- vector(mode = "list", length = length(snapshots))
  nPoints <- 0

  complete_owin <- spatstat.geom::owin(c(0,l), c(0,l))
  for (t_ind in 1:(length(snapshots)-1)) {
    complete_owin <- spatstat.geom::union.owin(complete_owin, spatstat.geom::owin(xrange = c(spacing*t_ind*l,(spacing*t_ind+1)*l), yrange = c(spacing*t_ind*l,(spacing*t_ind+1)*l)))
  }

  combined.points <- points_input[[1]][, 1:2]*l
  for (t_ind in 1:(length(snapshots)-1)) {
    t <- snapshots[t_ind]
    if(!is.numeric(points_input[[t]])){
      combined.points <- rbind(combined.points, points_input[[t]][, 1:2]*l + spacing*t_ind*l)
      nPoints <- nPoints + length(points_input[[t]]$x)
    }
  }
  l_hat <- nPoints/length(snapshots)
  combined.points <- spatstat.geom::ppp(x = combined.points$x, y = combined.points$y, window = complete_owin)

  return(spatstat.explore::Kest(combined.points, r = seq(from = 0, to = rMax, by = dr)))
}


#' Function for estimating the K-function given a list of point patterns
#' The point patterns are assumed to come from a time series of point patterns
#'
#' @param points_input List of point patterns used for estimation
#' @param l Side length of observation window
#' @param spacing How far apart are observations window in order to be independent?
#' @param rMax Max radius for which estimation is performed
#' @param dr How fine grid for estimation?
#' @param timescale How many indices apart do snapshots need to be in order to be independent?
#' @param divisor See documentation from the spatstat-function pcf.ppp
#'
#' @export
pcf_est.unions <- function(points_input, l, spacing, rMax, dr, timescale, divisor = "r"){
  snapshots <- seq(from = 1, to = length(points_input), by = ceiling(timescale))
  points.ppp <- vector(mode = "list", length = length(snapshots))
  nPoints <- 0

  complete_owin <- spatstat.geom::owin(c(0,l), c(0,l))
  for (t_ind in 1:(length(snapshots)-1)) {
    complete_owin <- spatstat.geom::union.owin(complete_owin,
                                               spatstat.geom::owin(xrange = c(spacing*t_ind*l,(spacing*t_ind+1)*l), yrange = c(spacing*t_ind*l,(spacing*t_ind+1)*l)))
  }

  combined.points <- points_input[[1]][, 1:2]*l
  for (t_ind in 1:(length(snapshots)-1)) {
    t <- snapshots[t_ind]
    if(!is.numeric(points_input[[t]])){
      combined.points <- rbind(combined.points, points_input[[t]][, 1:2]*l + spacing*t_ind*l)
      nPoints <- nPoints + length(points_input[[t]]$x)
    }
  }
  l_hat <- nPoints/length(snapshots)
  combined.points <- spatstat.geom::ppp(x = combined.points$x, y = combined.points$y, window = complete_owin)

  if(is.null(bw)){
    bw <- spatstat.explore::bw.pcf(combined.points, rmax = rMax)
  }
  return(spatstat.explore::pcf.ppp(combined.points, r = seq(from = 0, to = rMax, by = dr), bw = bw, divisor = divisor))
}
