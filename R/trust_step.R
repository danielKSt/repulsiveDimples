#' Do one iteration within trust region
#' @description
#' Iterate from params_sim using importance sampling on the realizations patternSim, simulated from params_sim.
#' The iteration is restricted to the ball centered on the parameters used for simulation,
#' with radius equal to the trust region radius.
#'
#' @param x_0 Parameter values used to simulate patternSim
#' @param delta Trust region radius
#' @param trust_function Trust function
#'
#' @export
trust_step <- function(x_0, delta, trust_function){
  n <- length(x_0)
  if(n == 1){
    x_k <- trust_line_steps(p_k = c(1), x_k = x_0, trust_function = trust_function,
                            subsection_count = 13, delta = delta, x_0 = x_0)
    return(list(x_star = x_k,
                f_pred = trust_function(x_k)))
  }
  p <- matrix(data = 0, nrow = n, ncol = n)
  diag(p) <- 1
  k <- 1
  x_k <- trust_line_steps(x_k = x_0, p_k = p[, 1], x_0 = x_0, delta = delta, trust_function = trust_function)
  hit_boundary <- (sqrt(sum(abs(x_k-x_0)^2)) > 0.99*delta)
  rotated <- FALSE
  while(!hit_boundary && k < 2*n){
    z_j <- matrix(data = 0, nrow = n, ncol = n + 1)
    z_j[, 1] <- x_k
    for(j in 1:n){
      z_j[, j+1] <- trust_line_steps(p_k = p[, j], x_k = z_j[, j], x_0 = x_0,
                                     delta = delta, trust_function = trust_function)
    }
    for(j in 1:(n-1)){
      p[, j] <- p[, j+1]
    }
    p[, n] <- z_j[, n + 1] - z_j[, 1]
    if(det(p)==0){
      if(rotated){
        break
      }
      p <- matrix(data = 0, nrow = n, ncol = n)
      diag(p) <- 1
      p[2, 1] <- -1
      p[1, 2] <- 1
      rotated <- TRUE
    }
    x_k <- trust_line_steps(x_k = z_j[, n + 1], p_k = p[, n], x_0 = x_0,
                            delta = delta, trust_function = trust_function)
    k <- k + 1
    hit_boundary <- (sqrt(sum(abs(x_k-x_0)^2)) > 0.99*delta)
  }

  return(list(x_star = x_k,
              f_pred = trust_function(x_k)))
}

#' Line minimization using interpolation with maximum reach
#' @description
#' Line optimization for the trust steps
#' @param p_k Search direction
#' @param x_k Search starting point
#' @param x_0 Center of trust region
#' @param delta Radius of trust region
#' @param subsection_count Number of bisections
#' @param trust_function Trust function
#' @param iterations How many times to subdivide the line
#'
#' @export
trust_line_steps <- function(p_k, x_k, x_0, delta, trust_function, subsection_count = 11, iterations = 4){
  alpha_range <- find_alpha_range(x_k, x_0, p_k, delta)
  alpha_k <- seq(from = alpha_range[1], to = alpha_range[2], length.out = subsection_count)
  z_k <- matrix(data = NA, nrow = length(x_k), ncol = subsection_count)
  for(i in 1:length(x_k)){
    z_k[i, ] <- x_k[i] + p_k[i] * alpha_k
  }
  f_vals <- apply(X = z_k, MARGIN = 2, FUN = trust_function)
  j_opt <- which.min(f_vals)
  for(j in 1:iterations){
    if(j_opt < 3){
      j_opt <- 3
    } else if(j_opt > subsection_count - 2){
      j_opt <- subsection_count - 2
    }
    alpha_k <- seq(from = alpha_k[j_opt-2], to = alpha_k[j_opt+2], length.out = subsection_count)
    for(i in 1:length(x_k)){
      z_k[i, ] <- x_k[i] + p_k[i] * alpha_k
    }
    f_vals <- apply(X = z_k, MARGIN = 2, FUN = trust_function)
    j_opt <- which.min(f_vals)
  }
  return(z_k[ ,j_opt])
}

#' Function to determine allowable range for alpha in line search
#' @description
#' We want to do a line search satisfying ||x_k+alpha*p_k - x_0||_2 <= delta
#'
#' This function solves the simple quadratic equation needed to find the range of
#' alpha values that satisfies thiat requirement.
#'
#' @param x_k Previous iterate, origin of line search
#' @param x_0 Midpoint of trust region
#' @param p_k Direction of line search
#' @param delta Radius of trust region
#' @export
find_alpha_range <- function(x_k, x_0, p_k, delta){
  a <- t(p_k)%*%p_k
  b <- 2*t(p_k)%*%(x_k-x_0)
  const <- t(x_k-x_0)%*%(x_k-x_0) - delta
  determ <- sqrt(b^2-4*a*const)
  up_lim <- (determ - b)/(2*a)
  low_lim <- (-1)*(determ + b)/(2*a)
  return(c(min(c(low_lim, up_lim)), max(c(low_lim, up_lim))))
}
