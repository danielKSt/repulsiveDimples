
#' Matern thinning
#'
#' @description
#' Perform a Matern type I, II or III thinning on a point pattern
#'
#' @details
#' Every point gets a uniform age mark, and in a close pair the older point has the
#' larger mark.
#' \itemize{
#'   \item Type I removes every point that has another point within `repulsionRange`.
#'   \item Type II removes every point that has an older point within `repulsionRange`,
#'   whether or not that older point is itself removed.
#'   \item Type III visits the points from oldest to youngest and removes a point only
#'   if an older point within `repulsionRange` has been retained.
#' }
#'
#' @param initialPattern DataFrame with the unthinned point pattern
#' @param repulsionRange Range of hard-core repulsion
#' @param xrange vector with max and min of x-axis in observation window
#' @param yrange vector with max and min of y-axis in observation window
#' @param thinningType Type of Matern thinning to perform, 1, 2 or 3.
#'
#' @export
matern.thinning <- function(initialPattern, repulsionRange, xrange, yrange, thinningType = 2){
  if(length(thinningType) != 1 || !(thinningType %in% 1:3)){
    stop("thinningType must be 1, 2 or 3")
  }
  n <- nrow(initialPattern)
  ageMark <- stats::runif(n = n)
  initialPattern$age <- ageMark
  df <- initialPattern[order(initialPattern$age), ]

  removePoint <- rep(0, n)
  if(n >= 2){
    pad <- max(repulsionRange, 1e-6)
    win <- spatstat.geom::owin(range(df$x) + c(-pad, pad),
                                range(df$y) + c(-pad, pad))
    pp <- spatstat.geom::ppp(df$x, df$y, window = win, check = FALSE)
    # Each close pair is listed once with i < j, and since df is sorted by age,
    # j is always the older point of the pair.
    pairs <- spatstat.geom::closepairs(pp, rmax = repulsionRange, twice = FALSE, what = "indices")
    if(thinningType == 1){
      removePoint[unique(c(pairs$i, pairs$j))] <- 1
    } else if(thinningType == 2){
      removePoint[unique(pairs$i)] <- 1
    } else {
      # Visiting from oldest to youngest means every older neighbour of a point has
      # already been settled by the time the point itself is reached.
      olderNeighbours <- split(pairs$j, factor(pairs$i, levels = seq_len(n)))
      for(k in sort(unique(pairs$i), decreasing = TRUE)){
        if(any(removePoint[olderNeighbours[[k]]] == 0)){
          removePoint[k] <- 1
        }
      }
    }
  }

  res <- df[which(removePoint==0), c(1,2)]
  res <- res[which(res$x > xrange[1]), ]
  res <- res[which(res$x < xrange[2]), ]
  res <- res[which(res$y > yrange[1]), ]
  res <- res[which(res$y < yrange[2]), ]
  return(res)
}

#' Thomas process with Matern thinning
#'
#' @description
#' Simulates a Thomas process and applies a Matern thinning of type `thinningType` to it
#'
#' @details
#' The Thomas process is simulated on the window enlarged by `repulsionRange` on every
#' side, so that points near the edge are thinned against the points just outside it.
#' This is exact for types I and II, where whether a point is retained depends only on
#' the points within `repulsionRange` of it. It is not exact for type III, where that
#' dependence reaches further through chains of close points.
#'
#' @param kappa See rThomas in spatstat
#' @param scale See rThomas in spatstat
#' @param mu See rThomas in spatstat
#' @param repulsionRange Range of hard-core repulsion
#' @param xlims xlim
#' @param ylims ylim
#' @param saveparents Logical value indicating whether to save the locations of the parent points as an attribute.
#' @param thinningType Type of Matern thinning, 1, 2 or 3, see \code{\link{matern.thinning}}.
#'
#' @export
rThomas_matern_thinned <- function(kappa, scale, mu, repulsionRange, xlims, ylims, saveparents = FALSE,
                                   thinningType = 2){
  xlims_un <- c(xlims[1] - repulsionRange, xlims[2] + repulsionRange)
  ylims_un <- c(ylims[1] - repulsionRange, ylims[2] + repulsionRange)
  # I've set algorithm = 'naive' due to some issues with the default for large simulation windows.
  unthinned <- spatstat.random::rThomas(kappa = kappa, scale = scale, mu = mu,
                                        win = spatstat.geom::owin(xlims_un, ylims_un),
                                        algorithm = "naive",
                                        saveparents = saveparents)

  if(saveparents){
    parents <- attr(unthinned, "parents")
    B_area <- spatstat.geom::area(parents$window)
    if(parents$n == 0){
      parents <- data.frame()
    } else {
      parents <- data.frame(parents)
    }
  }

  if(unthinned$n == 0){
    if(saveparents){ unthinned <- data.frame() }
    thinned <- data.frame()
  } else if(unthinned$n == 1){
    unthinned <- data.frame(x = unthinned$x, y = unthinned$y)
    pointInReducedBox <- (unthinned$x[1] > xlims[1]) && (unthinned$x[1] < xlims[2]) && (unthinned$y[1] > ylims[1]) && (unthinned$y[1] < ylims[2])
    if(pointInReducedBox){
      thinned <- unthinned
    } else {
      thinned <- data.frame()
    }
  } else {
    unthinned <- data.frame(x = unthinned$x, y = unthinned$y)
    if(repulsionRange == 0){
      thinned <- unthinned
    } else {
      thinned <- matern.thinning(initialPattern = unthinned, repulsionRange = repulsionRange,
                                 xrange = xlims, yrange = ylims, thinningType = thinningType)
    }
  }

  # Return the simulated result
  if(saveparents){
    return(list(parent = parents, daughter = unthinned, thinned = thinned,
                xlim_unthinned = xlims_un, ylim_unthinned = ylims_un, B_area = B_area,
                xlim_thinned = xlims, ylim_thinned = ylims))
  } else {
    return(thinned)
  }
}


#' Variance Gamma with Matern thinning
#'
#' @description
#' Simulates a Variance Gamma SNCP and applies a Matern thinning of type `thinningType` to it
#'
#'
#' @param kappa See rVarGamma in spatstat
#' @param scale See rVarGamma in spatstat
#' @param mu See rVarGamma in spatstat
#' @param nu See rVarGamma in spatstat
#' @param repulsionRange Range of hard-core repulsion
#' @param win Simulation window to be used
#' @param thinningType Type of Matern thinning, 1, 2 or 3, see \code{\link{matern.thinning}}.
#'
#' @export
rVarGamma_matern_thinned <- function(kappa, scale, mu, nu, repulsionRange, win, thinningType = 2){
  # I've set algorithm = 'naive' due to some issues with the default for large simulation windows.
  unthinned <- spatstat.random::rVarGamma(kappa = kappa, scale = scale, mu = mu, nu = nu,
                                          win = win, algorithm = "naive")
  if(unthinned$n > 1){
    unthinned <- data.frame(x = unthinned$x, y = unthinned$y)
    thinned <- matern.thinning(initialPattern = unthinned, repulsionRange = repulsionRange,
                               xrange = win$xrange, yrange = win$yrange, thinningType = thinningType)
    return(spatstat.geom::ppp(x = thinned$x, y = thinned$y, window = win))
  } else {
    return(unthinned)
  }
}
