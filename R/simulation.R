
#' Matern type II thinning
#'
#' @description
#' Perform a Matern type II thinning on a point pattern
#'
#' @param initialPattern DataFrame with the unthinned point pattern
#' @param repulsionRange Range of hard-core repulsion
#' @param xrange vector with max and min of x-axis in observation window
#' @param yrange vector with max and min of y-axis in observation window
#'
#' @export
matern.thinning <- function(initialPattern, repulsionRange, xrange, yrange){
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
    pairs <- spatstat.geom::closepairs(pp, rmax = repulsionRange, twice = FALSE, what = "indices")
    removePoint[unique(pairs$i)] <- 1
  }

  res <- df[which(removePoint==0), c(1,2)]
  res <- res[which(res$x > xrange[1]), ]
  res <- res[which(res$x < xrange[2]), ]
  res <- res[which(res$y > yrange[1]), ]
  res <- res[which(res$y < yrange[2]), ]
  return(res)
}

#' Thomas process with Matern II thinning
#'
#' @description
#' Simulates a Variance Gamma SNCP and applies a Matern II thinning to it
#'
#'
#' @param kappa See rThomas in spatstat
#' @param scale See rThomas in spatstat
#' @param mu See rThomas in spatstat
#' @param repulsionRange Range of hard-core repulsion
#' @param xlims xlim
#' @param ylims ylim
#' @param saveparents Logical value indicating whether to save the locations of the parent points as an attribute.
#'
#' @export
rThomas_matern_thinned <- function(kappa, scale, mu, repulsionRange, xlims, ylims, saveparents = FALSE){
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
                                 xrange = xlims, yrange = ylims)
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


#' Variance Gamma with Matern II thinning
#'
#' @description
#' Simulates a Variance Gamma SNCP and applies a Matern II thinning to it
#'
#'
#' @param kappa See rVarGamma in spatstat
#' @param scale See rVarGamma in spatstat
#' @param mu See rVarGamma in spatstat
#' @param nu See rVarGamma in spatstat
#' @param repulsionRange Range of hard-core repulsion
#' @param win Simulation window to be used
#'
#' @export
rVarGamma_matern_thinned <- function(kappa, scale, mu, nu, repulsionRange, win){
  # I've set algorithm = 'naive' due to some issues with the default for large simulation windows.
  unthinned <- spatstat.random::rVarGamma(kappa = kappa, scale = scale, mu = mu, nu = nu,
                                          win = win, algorithm = "naive")
  if(unthinned$n > 1){
    unthinned <- data.frame(x = unthinned$x, y = unthinned$y)
    thinned <- matern.thinning(initialPattern = unthinned, repulsionRange = repulsionRange,
                               xrange = win$xrange, yrange = win$yrange)
    return(spatstat.geom::ppp(x = thinned$x, y = thinned$y, window = win))
  } else {
    return(unthinned)
  }
}
