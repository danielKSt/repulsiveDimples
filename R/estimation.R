
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
  nPoints <- 0

  for (t_ind in seq_along(snapshots)) {
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
#' @description
#' Translation-corrected estimate of Ripley's K-function from independent snapshots
#' observed in the same square window \eqn{W=[0,l]^2}. This is the union approach: the
#' snapshots are treated as one pattern on a union \eqn{W^*} of \eqn{m} disjoint, well
#' separated copies of \eqn{W}, and the result is what
#' \code{spatstat.explore::Kest(..., correction = "translate")} would give on that union,
#' \deqn{\hat K(r)=\frac{|W^*|^2}{n(n-1)}\sum_{t}\sum_{i\neq j}
#' \frac{1\{\|x_i-x_j\|\leq r\}}{|W^*\cap W^*_{x_i-x_j}|},}
#' with the inner sum over pairs within snapshot \eqn{t}, \eqn{n} the total number of
#' points and \eqn{|W^*\cap W^*_v|=m(l-|v_1|)(l-|v_2|)}. It is computed per snapshot
#' with \code{\link{estimate_K_lambda_baseline}}, the same translation-corrected sum the
#' simulated ensembles use, rather than by spatstat on the union window, whose set
#' covariance spatstat only approximates on a pixel grid for non-rectangular windows.
#'
#' @param points_input List of point patterns used for estimation. Each element is a
#' `data.frame` with columns `x` and `y` in observation window coordinates, that is in
#' `[0, l]`, or a bare number for a snapshot holding no points.
#' @param l Side length of observation window
#' @param spacing Unused. The snapshots used to be laid out `spacing*l` apart on an
#' explicit union window; the estimate does not depend on the layout as long as the
#' copies are far enough apart that no pair crosses between them. Kept so that existing
#' calls still work.
#' @param rMax Max radius for which estimation is performed
#' @param dr How fine grid for estimation?
#' @param r_vec Vector to specify which radii to estimate for
#' @param timescale How many indices apart do snapshots need to be in order to be independent?
#'
#' @return An `fv` object with columns `r`, `theo` (the Poisson value \eqn{\pi r^2}) and
#' `trans` (the estimate), so that it can be passed directly to e.g.
#' \code{spatstat.model::thomas.estK}.
#'
#' @export
K_est.unions <- function(points_input, l, spacing, rMax = 5, dr = 0.05, r_vec = NULL, timescale){
  if(is.null(r_vec)){
    r_vec <- seq(from = 0, to = rMax, by = dr)
  }
  snapshots <- seq(from = 1, to = length(points_input), by = ceiling(timescale))

  n_total <- 0
  pair_sum <- numeric(length(r_vec))
  for (t in snapshots) {
    current <- points_input[[t]]
    if(is.numeric(current)){
      next
    }
    n_total <- n_total + nrow(current)
    pair_sum <- pair_sum + estimate_K_lambda_baseline(
      pattern = list(thinned = current[, 1:2], xlim_thinned = c(0, l), ylim_thinned = c(0, l)),
      r_vec = r_vec)$trans
  }
  if(n_total < 2){
    stop("K_est.unions needs at least two points in total.")
  }

  # pair_sum/m is the sum of 1/|W* intersect W*_v| over all ordered pairs, and the
  # squared intensity is n(n-1)/|W*|^2, as in Kest.
  m <- length(snapshots)
  area <- m*l^2
  K_trans <- (pair_sum/m)*area^2/(n_total*(n_total - 1))

  return(spatstat.explore::fv(data.frame(r = r_vec, theo = pi*r_vec^2, trans = K_trans),
                              argu = "r", ylab = quote(K(r)), valu = "trans", fmla = . ~ r,
                              alim = range(r_vec),
                              labl = c("r", "%s[pois](r)", "hat(%s)[trans](r)"),
                              desc = c("distance argument r", "theoretical Poisson %s",
                                       "translation-corrected estimate of %s"),
                              fname = "K"))
}


#' Function for estimating the pair correlation function given a list of point patterns
#' The point patterns are assumed to come from a time series of point patterns
#'
#' @param points_input List of point patterns used for estimation. Each element is a
#' `data.frame` with columns `x` and `y` in observation window coordinates, that is in
#' `[0, l]`, or a bare number for a snapshot holding no points.
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

  complete_owin <- spatstat.geom::owin(c(0,l), c(0,l))
  for (t_ind in seq_len(length(snapshots) - 1)) {
    complete_owin <- spatstat.geom::union.owin(complete_owin,
                                               spatstat.geom::owin(xrange = c(spacing*t_ind*l,(spacing*t_ind+1)*l), yrange = c(spacing*t_ind*l,(spacing*t_ind+1)*l)))
  }

  first_found <- FALSE
  t_first <- 0
  while((!first_found) && (t_first < length(snapshots))){
    t_first <- t_first + 1
    first_found <- !is.numeric(points_input[[snapshots[t_first]]])
  }

  combined.points <- points_input[[snapshots[t_first]]][, 1:2]
  for (t_ind in seq_len(length(snapshots) - t_first) + t_first) {
    t <- snapshots[t_ind]
    if(!is.numeric(points_input[[t]])){
      combined.points <- rbind(combined.points, points_input[[t]][, 1:2] + spacing*(t_ind-1)*l)
    }
  }
  combined.points <- spatstat.geom::ppp(x = combined.points$x, y = combined.points$y, window = complete_owin)

  if(is.null(bw)){
    bw <- spatstat.explore::bw.pcf(combined.points, rmax = rMax)
  }
  return(spatstat.explore::pcf.ppp(combined.points, r = seq(from = 0, to = rMax, by = dr), bw = bw, divisor = divisor, zerocor = zerocor))
}


#' Minimum inter-point distance for list of points
#'
#' @description
#' Smallest distance between any two points of the same snapshot, taken over all
#' snapshots. Used to bound the hard-core repulsion range of the Matern II thinning.
#'
#' The smallest pairwise distance within a pattern is the smallest of its
#' nearest-neighbour distances, so each snapshot is handled by one
#' \code{spatstat.geom::nndist} call rather than by an explicit loop over all
#' \eqn{O(n^2)} pairs.
#'
#' @param points_input List of point patterns used for estimation. Each element is
#' either the full output of \code{\link{rThomas_matern_thinned}} (the retained points
#' are then taken from its `thinned` element) or a `data.frame` with columns `x` and
#' `y`. Snapshots stored as a bare number, the empty-pattern sentinel used by
#' \code{\link{K_est.unions}}, are skipped.
#' @param l Side length of observation window. Also the value returned when no snapshot
#' holds two or more points.
#' @param timescale How many indices apart do snapshots need to be in order to be independent?
#'
#' @return The smallest inter-point distance found, or `l` if that is smaller.
#'
#' @export
find_min_dist <- function(points_input, timescale = 1, l){
  res <- l
  snapshots <- seq(from = 1, to = length(points_input), by = ceiling(timescale))
  for (t in snapshots) {
    current <- points_input[[t]]
    if(is.numeric(current)){
      next
    }
    if(!is.null(current$thinned)){
      current <- current$thinned
    }
    if(length(current$x) < 2){
      next
    }
    res <- min(c(res, spatstat.geom::nndist(current$x, current$y)))
  }
  return(res)
}


