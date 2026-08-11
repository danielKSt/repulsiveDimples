
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
#' @param r_vec Vector to specify which radii to estimate for
#' @param timescale How many indices apart do snapshots need to be in order to be independent?
#'
#' @export
K_est.unions <- function(points_input, l, spacing, rMax = 5, dr = 0.05, r_vec = NULL, timescale){
  snapshots <- seq(from = 1, to = length(points_input), by = ceiling(timescale))
  points.ppp <- vector(mode = "list", length = length(snapshots))
  nPoints <- 0

  complete_owin <- spatstat.geom::owin(c(0,l), c(0,l))
  for (t_ind in 1:(length(snapshots)-1)) {
    complete_owin <- spatstat.geom::union.owin(complete_owin, spatstat.geom::owin(xrange = c(spacing*t_ind*l,(spacing*t_ind+1)*l), yrange = c(spacing*t_ind*l,(spacing*t_ind+1)*l)))
  }

  first_found <- FALSE
  t_first <- 0
  while((!first_found) && (t_first < length(snapshots))){
    t_first <- t_first + 1
    first_found <- !is.numeric(points_input[[snapshots[t_first]]])
  }

  combined.points <- points_input[[snapshots[t_first]]][, 1:2]*l
  for (t_ind in (t_first+1):length(snapshots)) {
    t <- snapshots[t_ind]
    if(!is.numeric(points_input[[t]])){
      combined.points <- rbind(combined.points, points_input[[t]][, 1:2]*l + spacing*(t_ind-1)*l)
      nPoints <- nPoints + length(points_input[[t]]$x)
    }
  }
  l_hat <- nPoints/length(snapshots)
  combined.points <- spatstat.geom::ppp(x = combined.points$x, y = combined.points$y, window = complete_owin)

  if(is.null(r_vec)){
    return(spatstat.explore::Kest(combined.points, r = seq(from = 0, to = rMax, by = dr), correction = c("border")))
  } else {
    return(spatstat.explore::Kest(combined.points, r = r_vec, correction = c("border")))
  }
}


#' Function for estimating the K-function given a list of point patterns
#' The point patterns are assumed to come from a time series of point patterns
#'
#' @param points_input List of point patterns used for estimation
#' @param bw Estimation bandwidth
#' @param l Side length of observation window
#' @param spacing How far apart are observations window in order to be independent?
#' @param rMax Max radius for which estimation is performed
#' @param dr How fine grid for estimation?
#' @param timescale How many indices apart do snapshots need to be in order to be independent?
#' @param divisor See documentation from the spatstat-function pcf.ppp
#' @param zerocor See pcf.ppp documentation from spatstat
#'
#' @export
pcf_est.unions <- function(points_input, bw = NULL, l, spacing, rMax, dr, timescale, divisor = "r", zerocor = NULL){
  snapshots <- seq(from = 1, to = length(points_input), by = ceiling(timescale))
  points.ppp <- vector(mode = "list", length = length(snapshots))
  nPoints <- 0

  complete_owin <- spatstat.geom::owin(c(0,l), c(0,l))
  for (t_ind in 1:(length(snapshots)-1)) {
    complete_owin <- spatstat.geom::union.owin(complete_owin,
                                               spatstat.geom::owin(xrange = c(spacing*t_ind*l,(spacing*t_ind+1)*l), yrange = c(spacing*t_ind*l,(spacing*t_ind+1)*l)))
  }

  first_found <- FALSE
  t_first <- 0
  while((!first_found) && (t_first < length(snapshots))){
    t_first <- t_first + 1
    first_found <- !is.numeric(points_input[[snapshots[t_first]]])
  }

  combined.points <- points_input[[snapshots[t_first]]][, 1:2]*l
  for (t_ind in (t_first+1):length(snapshots)) {
    t <- snapshots[t_ind]
    if(!is.numeric(points_input[[t]])){
      combined.points <- rbind(combined.points, points_input[[t]][, 1:2]*l + spacing*(t_ind-1)*l)
      nPoints <- nPoints + length(points_input[[t]]$x)
    }
  }
  l_hat <- nPoints/length(snapshots)
  combined.points <- spatstat.geom::ppp(x = combined.points$x, y = combined.points$y, window = complete_owin)

  if(is.null(bw)){
    bw <- spatstat.explore::bw.pcf(combined.points, rmax = rMax)
  }
  return(spatstat.explore::pcf.ppp(combined.points, r = seq(from = 0, to = rMax, by = dr), bw = bw, divisor = divisor, zerocor = zerocor))
}


#' Minimum inter-point distance for list of points
#'
#' @param points_input List of point patterns used for estimation
#' @param l Side length of observation window
#' @param timescale How many indices apart do snapshots need to be in order to be independent?
#'
#' @export
find_min_dist <- function(points_input, timescale = 1, l){
  res <- l
  snapshots <- seq(from = 1, to = length(points_input), by = ceiling(timescale))
  for (t in 1:length(snapshots)) {
    current <- points_input[[snapshots[t]]]$thinned
    if(nrow(current) > 2){
      for(i in 1:(length(current$x)-1)){
        for (j in (i+1):length(current$x)) {
          dist <- sqrt((current$x[i]-current$x[j])^2+(current$y[i]-current$y[j])^2)
          res <- min(c(dist,res))
        }
      }
    }
  }
  return(res)
}


