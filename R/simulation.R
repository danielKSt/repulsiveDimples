
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
  unthinned <- spatstat.random::rVarGamma(kappa = kappa, scale = scale, mu = mu, nu = nu, win = win, algorithm = "naive")
  unthinned <- data.frame(x = unthinned$x, y = unthinned$y)
  thinned <- matern.thinning(initialPattern = unthinned, repulsionRange = repulsionRange, xrange = win$xrange, yrange = win$yrange)
  return(spatstat.geom::ppp(x = thinned$x, y = thinned$y, window = win))
}
