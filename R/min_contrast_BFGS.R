

#' Complete estimation procedure using importance sampling + simulation + minimum contrast
#' @description
#' Minimum contrast estimation using trust region approach with simulation based estimation
#' of gradient and SR1 for Hessian approximation.
#' This method also uses importance sampling...
#'
#'
#' @param tol Convergence tolerance
#' @param max.iter Maximum iterations
#' @param params Initial guess for parameters to estimate
#' @param par_free_index Indices of parameters to estimate
#' @param par_fixed Values of fixed parameters
#' @param epsilon Step for derivative calculation
#' @param imp_samp_tolerance Tolerance for difference in parameter value before resampling patternSim
#' @param delta_hat trust radius
#' @param simulate_point_pattern Function for simulation of point pattern
#'
#' @export
minimum_contrast_estimation <- function(params, par_free_index, par_fixed, epsilon, imp_samp_tolerance, delta_hat,
                                        tol = 10^-8, max.iter = 10000, simulate_point_pattern){
  # Initialize output:
  params_sequence <- matrix(data = 0, nrow = max.iter, ncol = length(params))
  f_vals <- rep(0, max.iter)
  gn_vals <- matrix(data = 0, nrow = max.iter, ncol = length(params))

  # Initaial step
  patternSim <- simulate_point_pattern(params, par_free_index, par_fixed)
  f_vals[1] <- contrast_function_importance_sampling(patternSim, params, par_free_index, par_fixed)
  gn_vals[1, ] <- sapply(X = c(1:length(params)), FUN = contrast_gradient_importance_sampling,
                         par_free_index = par_free_index,
                         par_fixed = par_fixed,
                         params = params,
                         patternSim = patternSim,
                         epsilon = epsilon,
                         f_val = f_vals[1])
  H <- matrix(data = 1, nrow = 2, ncol = 2)
  #H <- contrast_H_init_importance_sampling(patternSim, params, par_free_index, par_fixed, epsilon)
  #params_sequence[1, ] <- trust_region_update(...)
  params_sequence[1, ] <- c(1,1,1)
  nSteps <- 1
  converged <- (gn_vals[1, ]%*%gn_vals[1, ] < tol)
  while((nSteps < max.iter) && !converged){
    patternSim <- simulate_point_pattern(params, par_free_index, par_fixed)
    f_vals[nSteps+1] <- contrast_function_importance_sampling(patternSim, params, par_free_index, par_fixed)
    gn_vals[nSteps+1, ] <- sapply(X = c(1:length(params)), FUN = contrast_gradient_importance_sampling,
                                par_free_index = par_free_index,
                                par_fixed = par_fixed,
                                params = params,
                                patternSim = patternSim,
                                epsilon = epsilon,
                                f_val = f_vals[nSteps+1])

    nSteps <- nSteps + 1
  }

}

#' Calculate partial derivative for the k'th parameter
#' @description
#' Calculate partial derivative for the k'th parameter
#'
#' @param k Index of direction to estimate derivative in
#' @param patternSim simulated point pattern to use
#' @param epsilon Step for derivative calculation
#' @param params Parameter values at which to aproximate the partial derivative
#' @param par_free_index Indices of parameters to estimate
#' @param par_fixed Values of fixed parameters
#' @param f_val Value of function at evaluation location
#'
#' @export
contrast_gradient_importance_sampling <- function(params, par_free_index, par_fixed, k, patternSim, epsilon, f_val){
  step_params <- rep(0, length(params))
  step_params[k] <- epsilon
  f_fwd <- contrast_function_importance_sampling(patternSim, params, par_free_index, par_fixed)
  return((f_fwd - f_val)/epsilon)
}

#' Perform a step of the trust region update
#' @description
#' A short description...
#'
#' @export
trust_region_update <- function(){
  return(0)
}





























