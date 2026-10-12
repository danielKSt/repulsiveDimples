
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
#' The parents are simulated with `nonempty = FALSE`, so the saved parents are a Poisson
#' process of intensity `kappa`, parents without daughters included. The importance
#' sampling weights assume this. spatstat's default, `nonempty = TRUE`, only keeps parents
#' with at least one daughter, which leaves the daughters' distribution unchanged but
#' makes the raw importance sampling weights biased by a constant factor.
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
#' @return If `saveparents` is `FALSE`, a data frame with the thinned points. Otherwise the
#' list described in \code{\link{matern_thin_cluster_pattern}}.
#'
#' @export
rThomas_matern_thinned <- function(kappa, scale, mu, repulsionRange, xlims, ylims, saveparents = FALSE,
                                   thinningType = 2){
  xlims_un <- c(xlims[1] - repulsionRange, xlims[2] + repulsionRange)
  ylims_un <- c(ylims[1] - repulsionRange, ylims[2] + repulsionRange)
  # I've set algorithm = 'naive' due to some issues with the default for large simulation windows.
  unthinned <- spatstat.random::rThomas(kappa = kappa, scale = scale, mu = mu,
                                        win = spatstat.geom::owin(xlims_un, ylims_un),
                                        algorithm = "naive", nonempty = FALSE,
                                        saveparents = saveparents)

  return(matern_thin_cluster_pattern(unthinned = unthinned, repulsionRange = repulsionRange,
                                     xlims = xlims, ylims = ylims, saveparents = saveparents,
                                     thinningType = thinningType))
}


#' Variance Gamma with Matern thinning
#'
#' @description
#' Simulates a Variance Gamma SNCP and applies a Matern thinning of type `thinningType` to it
#'
#' @details
#' Simulated in the same way as \code{\link{rThomas_matern_thinned}}: on the window enlarged
#' by `repulsionRange`, and with `nonempty = FALSE`, see there for why.
#'
#' The cluster kernel is
#' \deqn{k(u) = \frac{(\|u\|/\sigma)^{\nu} K_{\nu}(\|u\|/\sigma)}{\pi 2^{\nu+1} \sigma^2 \Gamma(\nu+1)}}
#' with \eqn{\sigma} = `scale`, where \eqn{K_\nu} is the modified Bessel function of the
#' second kind.
#'
#' @param kappa See rVarGamma in spatstat
#' @param scale See rVarGamma in spatstat
#' @param mu See rVarGamma in spatstat
#' @param nu See rVarGamma in spatstat
#' @param repulsionRange Range of hard-core repulsion
#' @param xlims xlim
#' @param ylims ylim
#' @param saveparents Logical value indicating whether to save the locations of the parent points as an attribute.
#' @param thinningType Type of Matern thinning, 1, 2 or 3, see \code{\link{matern.thinning}}.
#'
#' @return If `saveparents` is `FALSE`, a data frame with the thinned points. Otherwise the
#' list described in \code{\link{matern_thin_cluster_pattern}}.
#'
#' @export
rVarGamma_matern_thinned <- function(kappa, scale, mu, nu, repulsionRange, xlims, ylims, saveparents = FALSE,
                                     thinningType = 2){
  xlims_un <- c(xlims[1] - repulsionRange, xlims[2] + repulsionRange)
  ylims_un <- c(ylims[1] - repulsionRange, ylims[2] + repulsionRange)
  # I've set algorithm = 'naive' due to some issues with the default for large simulation windows.
  unthinned <- spatstat.random::rVarGamma(kappa = kappa, scale = scale, mu = mu, nu = nu,
                                          win = spatstat.geom::owin(xlims_un, ylims_un),
                                          algorithm = "naive", nonempty = FALSE,
                                          saveparents = saveparents)

  return(matern_thin_cluster_pattern(unthinned = unthinned, repulsionRange = repulsionRange,
                                     xlims = xlims, ylims = ylims, saveparents = saveparents,
                                     thinningType = thinningType))
}


#' Matern thinning of a simulated cluster process
#'
#' @description
#' The step shared by \code{\link{rThomas_matern_thinned}} and
#' \code{\link{rVarGamma_matern_thinned}}: applies the Matern thinning to a cluster process
#' simulated on the window enlarged by `repulsionRange`, keeps the thinned points inside
#' `xlims` x `ylims`, and collects the parents and unthinned daughters if they were saved.
#'
#' @param unthinned The unthinned cluster process, a \code{ppp} simulated on the window
#' enlarged by `repulsionRange` on every side, with the parents as an attribute if
#' `saveparents` is `TRUE`.
#' @param repulsionRange Range of hard-core repulsion
#' @param xlims xlim of the observation window
#' @param ylims ylim of the observation window
#' @param saveparents Logical value indicating whether the parents were saved.
#' @param thinningType Type of Matern thinning, 1, 2 or 3, see \code{\link{matern.thinning}}.
#'
#' @return If `saveparents` is `FALSE`, a data frame with the thinned points. Otherwise a
#' list with `parent` (the parents in the enlarged parent window \eqn{B}), `daughter` (the
#' unthinned daughters in the enlarged window), `thinned`, `xlim_unthinned` and
#' `ylim_unthinned` (the enlarged window), `B_area`, and `xlim_thinned` and `ylim_thinned`.
matern_thin_cluster_pattern <- function(unthinned, repulsionRange, xlims, ylims, saveparents, thinningType){
  xlims_un <- c(xlims[1] - repulsionRange, xlims[2] + repulsionRange)
  ylims_un <- c(ylims[1] - repulsionRange, ylims[2] + repulsionRange)

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
