
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
  # Bruk closePairs for å forbetre
  # Kutt ut ytre kant av vindu i sluttresultatet
  ageMark <- stats::runif(n = nrow(initialPattern))
  initialPattern$age <- ageMark
  df <- initialPattern[order(initialPattern$age), ]
  removePoint <- rep(0, nrow(df))
  for(i in 1:(nrow(df)-1)){
    min_dist = sqrt(xrange[2]^2+yrange[2]^2)
    for (j in (i+1):nrow(df)) {
      d <- sqrt((df$x[i]-df$x[j])^2+(df$y[i]-df$y[j])^2)
      min_dist <- min(c(d, min_dist))
      if(min_dist <= repulsionRange){
        removePoint[i] <- 1
        break
      }
    }
  }
  res <- df[which(removePoint==0), c(1,2)]
  return(res)
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
#' @param win Simulation window to be used
#' @param saveparents Logical value indicating whether to save the locations of the parent points as an attribute.
#'
#' @export
rThomas_matern_thinned <- function(kappa, scale, mu, repulsionRange, win, saveparents = FALSE){
  # I've set algorithm = 'naive' due to some issues with the default for large simulation windows.
  unthinned <- spatstat.random::rThomas(kappa = kappa, scale = scale, mu = mu,
                                        win = win, algorithm = "naive",
                                        saveparents = saveparents)

  expand_parent <- 4*scale
  B_area <- (win$xrange[2] - win$xrange[1]+2*expand_parent)*(win$yrange[2] - win$yrange[1]+2*expand_parent)
  if(unthinned$n > 1){
    if(saveparents){
      parents <- attr(unthinned, "parents")
      parents <- data.frame(parents)
    }
    unthinned <- data.frame(x = unthinned$x, y = unthinned$y)
    if(repulsionRange == 0){
      thinned <- spatstat.geom::ppp(x = unthinned$x, y = unthinned$y, window = win)
    } else {
      thinned <- matern.thinning(initialPattern = unthinned, repulsionRange = repulsionRange,
                                 xrange = win$xrange, yrange = win$yrange)
      thinned <- spatstat.geom::ppp(x = thinned$x, y = thinned$y, window = win)
    }
    if(saveparents){
      return(list(parent = parents, daughter = unthinned, thinned = thinned,
                  xlim = win$xrange, ylim = win$yrange, B_area = B_area))
    } else {
      return(thinned)
    }
  } else {
    if(saveparents){
      parents <- attr(unthinned, "parents")
      parents <- data.frame(parents)
      return(list(parent = parents, daughter = unthinned, thinned = unthinned,
                  xlim = win$xrange, ylim = win$yrange, B_area = B_area))
    } else {
      return(unthinned)
    }
  }
}
