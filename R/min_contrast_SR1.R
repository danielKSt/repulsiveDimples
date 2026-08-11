

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
#' @param epsilon Step for derivative calculation
#' @param delta_hat trust radius
#' @param K_hat Estimated K-function for data
#' @param rho_hat Estimated intensity for data
#' @param wq Weights for contrast function
#' @param normalized Normalize IS weights?
#' @param nSims Number of simulations
#' @param repRange Repulsion range
#' @param xlims Simulation window x limits
#' @param ylims Simulation window y limits
#' @param wq Weights for contrast function
#' @param eta Trust region parameter
#' @param simulation_threshold simulation_threshold
#' @param r Hyperparameter for SR1 update
#' @param printProgress Set to TRUE to get updates on progress while running
#'
#' @export
minimum_contrast_sr1 <- function(params, par_free_index, epsilon, delta_hat,
                                        tol = 10^-8, max.iter = 10000, normalized,
                                        rho_hat, K_hat, nSims, repRange, xlims, ylims, wq, eta,
                                        simulation_threshold, r, printProgress = FALSE){
  # Initialize output: ----
  nSteps <- 1
  params_sequence <- matrix(data = 0, nrow = max.iter, ncol = length(params))
  f_vals <- rep(0, max.iter)
  gn_vals <- matrix(data = 0, nrow = max.iter, ncol = length(par_free_index))

  # Initaial step ----
  if(printProgress){
    initTime <- Sys.time()
    print("Starting initial step")

    print("Simulating pattern: ")
    patternSim <- mcprogress::pmclapply(X = rep(exp(params[1]), nSims), FUN = rThomas_matern_thinned,
                            scale = exp(params[2]), mu = exp(params[3]),
                            repulsionRange = repRange, xlims = xlims, ylims = ylims, saveparents = TRUE)
    print("Estimating baselines: ")
    K_lambda_baseline <- mcprogress::pmclapply(patternSim, estimate_K_lambda_baseline, r_vec = K_hat$r)
  } else {
    patternSim <- parallel::mclapply(X = rep(exp(params[1]), nSims), FUN = rThomas_matern_thinned,
                                     scale = exp(params[2]), mu = exp(params[3]),
                                     repulsionRange = repRange, xlims = xlims, ylims = ylims, saveparents = TRUE)
    K_lambda_baseline <- parallel::mclapply(patternSim, estimate_K_lambda_baseline, r_vec = K_hat$r)
  }

  params_sim <- params
  rho_baseline <- sapply(patternSim, estimate_rho_baseline)

  evals <- evaluate_contrast_and_gradient(patternSim = patternSim, params = params, params_0 = params,
                                         rho_baseline = rho_baseline, K_lambda_baseline = K_lambda_baseline,
                                         rho_hat = rho_hat, K_hat = K_hat,
                                         par_free_index = par_free_index, epsilon = epsilon,
                                         normalized = normalized, wq = wq)

  f_vals[1] <- evals$f_val
  gn_vals[1, ] <- evals$gn_val[par_free_index]

  H <- contrast_hessian_importance_sampling(params = params, epsilon = epsilon,
                                            par_free_index = par_free_index, patternSim = patternSim,
                                            rho_baseline = rho_baseline,
                                            K_lambda_baseline = K_lambda_baseline,
                                            normalized = normalized,
                                            rho_hat = rho_hat, K_hat = K_hat,
                                            wq = wq, f_val = f_vals[1])

  trust_region_step <- quadratic_trust_region_update(evals$gn_val[par_free_index], H, trust_radius = delta_hat, tol = tol)

  params_new <- params
  params_new[par_free_index] <- params[par_free_index] + trust_region_step$params

  update_evaluation <- evaluate_improvement_SR1(f_old = f_vals[1], params_new = params_new, params_sim = params_sim,
                                            params_old = params, predicted_reduction = trust_region_step$predicted_red,
                                            patternSim = patternSim, K_hat = K_hat, rho_hat = rho_hat, wq = wq,
                                            rho_baseline = rho_baseline, K_lambda_baseline = K_lambda_baseline,
                                            normalized = normalized, delta_hat = delta_hat, eta = eta)

  if(update_evaluation$updateIterate){
    params_sequence[1, ] <- params_new
  } else {
    params_sequence[1, ] <- params
  }

  delta_hat <- update_evaluation$newDelta
  converged <- (gn_vals[1, ]%*%gn_vals[1, ] < tol)

  if(printProgress){
    print(paste("Initial iteration completed using: ", Sys.time() - initTime))
  }

  # Iteration loop ----
  while((nSteps < max.iter) && !converged){
    params <- params_sequence[nSteps, ]
    nSteps <- nSteps + 1
    if(printProgress && (nSteps %% 20 == 0)){
      endTime <- Sys.time()
      print(paste("20 steps completed in: ", endTime - initTime))
      initTime <- Sys.time()
      print(paste("Total steps completed: ", nSteps))
    }
    if(sqrt(sum((params - params_sim)^2)) > simulation_threshold){
      if(printProgress){
        initTime <- Sys.time()
        print(paste("Simulating new pattern at step", nSteps, ": "))
        patternSim <- mcprogress::pmclapply(X = rep(exp(params[1]), nSims), FUN = rThomas_matern_thinned,
                                            scale = exp(params[2]), mu = exp(params[3]),
                                            repulsionRange = repRange, xlims = xlims, ylims = ylims, saveparents = TRUE)
        print("Estimating baselines: ")
        K_lambda_baseline <- mcprogress::pmclapply(patternSim, estimate_K_lambda_baseline, r_vec = K_hat$r)
      } else {
        patternSim <- parallel::mclapply(X = rep(exp(params[1]), nSims), FUN = rThomas_matern_thinned,
                                         scale = exp(params[2]), mu = exp(params[3]),
                                         repulsionRange = repRange, xlims = xlims, ylims = ylims, saveparents = TRUE)
        K_lambda_baseline <- parallel::mclapply(patternSim, estimate_K_lambda_baseline, r_vec = K_hat$r)
      }

      rho_baseline <- sapply(patternSim, estimate_rho_baseline)
      params_sim <- params
    }

    if(update_evaluation$updateIterate){
      evals <- evaluate_contrast_and_gradient(patternSim = patternSim, params = params, params_0 = params_sim,
                                              rho_baseline = rho_baseline, K_lambda_baseline = K_lambda_baseline,
                                              rho_hat = rho_hat, K_hat = K_hat,
                                              par_free_index = par_free_index, epsilon = epsilon,
                                              normalized = normalized, wq = wq)

      f_vals[nSteps] <- evals$f_val
      gn_vals[nSteps, ] <- evals$gn_val[par_free_index]
    } else {
      f_vals[nSteps] <- f_vals[nSteps - 1]
      gn_vals[nSteps, ] <- gn_vals[nSteps - 1, ]
    }

    H <- SR1_update(Bk = H, sk = trust_region_step$params, r = r, yk = gn_vals[nSteps, ]- gn_vals[nSteps - 1, ])

    trust_region_step <- quadratic_trust_region_update(evals$gn_val[par_free_index], H, trust_radius = delta_hat, tol = tol)

    params_new <- params
    params_new[par_free_index] <- params[par_free_index] + trust_region_step$params

    update_evaluation <- evaluate_improvement_SR1(f_old = f_vals[nSteps], params_new = params_new, params_sim = params_sim,
                                              params_old = params, predicted_reduction = trust_region_step$predicted_red,
                                              patternSim = patternSim, K_hat = K_hat, rho_hat = rho_hat, wq = wq,
                                              rho_baseline = rho_baseline, K_lambda_baseline = K_lambda_baseline,
                                              normalized = normalized, delta_hat = delta_hat, eta = eta)

    if(update_evaluation$updateIterate){
      params_sequence[nSteps, ] <- params_new
    } else {
      params_sequence[nSteps, ] <- params_sequence[nSteps - 1, ]
    }

    delta_hat <- update_evaluation$newDelta
    converged <- (gn_vals[nSteps, ]%*%gn_vals[nSteps, ] < tol)
  }
  return(list(params = params_sequence,
              f_vals = f_vals,
              gradients = gn_vals))
}


#' Evaluate if there is sufficient improvement
#' @description
#' evaluate improvement
#'
#' @param f_old Value of function at evaluation location
#' @param params_new Parameter values at which to approximate the partial derivative
#' @param params_old Parameter values of the previous iteration
#' @param params_sim Parameter values used to simulate  patternSim
#' @param predicted_reduction Step for derivative calculation
#' @param patternSim simulated point pattern to use
#' @param rho_baseline rho_baseline
#' @param K_lambda_baseline K_lambda_baseline
#' @param normalized Normalize IS weights?
#' @param rho_hat Estimated intensity for data
#' @param K_hat Estimated K-function for data
#' @param wq Weights for contrast function
#' @param delta_hat Old trust region radius
#' @param eta Trust region parameter
#'
#' @export
evaluate_improvement_SR1 <- function(f_old, params_new, params_old, params_sim, predicted_reduction, patternSim, K_hat, rho_hat, wq,
                                 rho_baseline, K_lambda_baseline, normalized, delta_hat, eta){
  f_new <- contrast_is(patternSim = patternSim, params_0 = params_sim, params_new = params_new,
                       K_hat = K_hat, rho_hat = rho_hat, wq = wq, rho_baseline = rho_baseline,
                       K_lambda_baseline = K_lambda_baseline, normalized = normalized)$f_est
  updateIterate <- ((f_old - f_new)/(predicted_reduction) > eta)
  if((f_old - f_new)/(predicted_reduction) > 0.75){
    if(sqrt(sum((params_new - params_old)^2)) > 0.8){
      delta_hat <- 2*delta_hat
    }
  } else if((f_old - f_new)/(predicted_reduction) < 0.1){
    delta_hat <- 0.5*delta_hat
  }
  return(list(updateIterate = updateIterate, newDelta = delta_hat))
}

#' Estimate gradient at some parameter value using importance sampling
#' @description
#' Calculate partial derivative for the k'th parameter
#'
#' @param params Parameter values at which to aproximate the partial derivative
#' @param params_sim Parameter values used in simulating patternSim
#' @param par_free_index Indices of parameters to estimate
#' @param epsilon Step for derivative calculation
#' @param patternSim simulated point pattern to use
#' @param rho_baseline rho_baseline
#' @param K_lambda_baseline K_lambda_baseline
#' @param normalized Normalize IS weights?
#' @param rho_hat Estimated intensity for data
#' @param K_hat Estimated K-function for data
#' @param wq Weights for contrast function
#' @param f_val Value of function at evaluation location
#'
#' @export
contrast_gradient_importance_sampling <- function(params, params_sim, par_free_index, epsilon, patternSim,
                                                  rho_baseline, K_lambda_baseline, normalized,
                                                  rho_hat, K_hat, wq, f_val){
  gn <- rep(0, length(params))
  for (i in 1:length(par_free_index)) {
    step_params <- params
    step_params[par_free_index[i]] <- params[par_free_index[i]] + epsilon

    f_fwd <- contrast_is(patternSim = patternSim,
                         params_0 = params_sim,
                         params_new = step_params,
                         rho_hat = rho_hat,
                         K_hat = K_hat,
                         rho_baseline = rho_baseline,
                         K_lambda_baseline = K_lambda_baseline,
                         normalized = normalized,
                         wq = wq)$f_est

    gn[par_free_index[i]] <- (f_fwd - f_val)/epsilon
  }
  return(gn)
}


#' Estimate Hessian at some parameter value using importance sampling
#' @description
#' See equation 8.21 in Nocedal and Wright.
#'
#' @param params Parameter values at which to aproximate the partial derivative
#' @param par_free_index Indices of parameters to estimate
#' @param epsilon Step for derivative calculation
#' @param patternSim simulated point pattern to use
#' @param rho_baseline rho_baseline
#' @param K_lambda_baseline K_lambda_baseline
#' @param normalized Normalize IS weights?
#' @param rho_hat Estimated intensity for data
#' @param K_hat Estimated K-function for data
#' @param wq Weights for contrast function
#' @param f_val Value of function at evaluation location
#'
#' @export
contrast_hessian_importance_sampling <- function(params, par_free_index, epsilon, patternSim,
                                                 rho_baseline, K_lambda_baseline, normalized,
                                                 rho_hat, K_hat, wq, f_val){
  H <- matrix(data = 1, nrow = length(par_free_index), ncol = length(par_free_index))
  for (i in 1:length(par_free_index)) {
    for(j in 1:length(par_free_index)) {
      step_params <- params
      step_params[i] <- params[i] + epsilon
      f_fwd_i <- contrast_is(patternSim = patternSim,
                           params_0 = params,
                           params_new = step_params,
                           rho_hat = rho_hat,
                           K_hat = K_hat,
                           rho_baseline = rho_baseline,
                           K_lambda_baseline = K_lambda_baseline,
                           normalized = normalized,
                           wq = wq)$f_est

      step_params <- params
      step_params[j] <- params[j] + epsilon
      f_fwd_j <- contrast_is(patternSim = patternSim,
                             params_0 = params,
                             params_new = step_params,
                             rho_hat = rho_hat,
                             K_hat = K_hat,
                             rho_baseline = rho_baseline,
                             K_lambda_baseline = K_lambda_baseline,
                             normalized = normalized,
                             wq = wq)$f_est

      step_params[i] <- step_params[i] + epsilon
      f_fwd_ij <- contrast_is(patternSim = patternSim,
                             params_0 = params,
                             params_new = step_params,
                             rho_hat = rho_hat,
                             K_hat = K_hat,
                             rho_baseline = rho_baseline,
                             K_lambda_baseline = K_lambda_baseline,
                             normalized = normalized,
                             wq = wq)$f_est

      H[i,j] <- (f_fwd_ij - f_fwd_i - f_fwd_j + f_val)/epsilon
    }
  }
  return(H)
}


#' Estimate both the contrast function and gradient of the contrast function using importance sampling
#' @description
#' ...?
#'
#' @param patternSim Simulated pattern
#' @param params Parameters to use in estimation
#' @param params_0 Parameters used for simulation
#' @param rho_baseline rho_baseline
#' @param K_lambda_baseline K_lambda_baseline
#' @param rho_hat Estimated intensity for data
#' @param K_hat Estimated K-function for data
#' @param par_free_index Indices of parameters to estimate
#' @param epsilon Step for derivative calculation
#' @param normalized Normalize IS weights?
#' @param wq Weights for contrast function
#'
#' @export
evaluate_contrast_and_gradient <- function(patternSim, params, params_0, rho_baseline, K_lambda_baseline, rho_hat, K_hat,
                                           par_free_index, epsilon, normalized, wq){

  if(is.null(params_0)){
    rho_par <- mean(rho_baseline)
    K_par <- K_importance_sampling(w_is = rep(1, length(rho_baseline)), patternSim = NULL,
                                   K_lambda_baseline = K_lambda_baseline,
                                   rho_baseline = rho_baseline)
    f_val <- contrast_function(rho_par = rho_par, rho_hat = rho_hat, K_par = K_par, K_hat = K_hat)
  } else {
    f_val <- contrast_is(patternSim = patternSim, params_0 = params_0, params_new = params,
                         K_hat = K_hat, rho_hat = rho_hat, wq = wq, rho_baseline = rho_baseline,
                         K_lambda_baseline = K_lambda_baseline, normalized = normalized)$f_est
  }

  gn_val <- contrast_gradient_importance_sampling(params = params, params_sim = params_0, par_free_index = par_free_index,
                                                  epsilon = epsilon, patternSim = patternSim,
                                                  rho_baseline = rho_baseline, K_lambda_baseline = K_lambda_baseline,
                                                  normalized = normalized,
                                                  rho_hat = rho_hat, K_hat = K_hat,
                                                  wq = wq, f_val = f_val)

  return(list(f_val = f_val, gn_val = gn_val))
}


#' Solve the trust-region subproblem via Steihaug-Toint truncated CG
#'
#' Solves:   min  g'p + (1/2) p'Hp
#'           s.t. ||p|| <= trust_radius
#'
#' The Steihaug-Toint method extends conjugate gradients with two extra
#' early-exit conditions, making it well-suited for SR1 trust-region methods
#' where the Hessian approximation H may be indefinite.
#'
#' Each CG iteration ends in one of three ways:
#'   (1) Negative/zero curvature (d'Hd <= 0): step to the trust-region
#'       boundary along d, choosing the intersection point with lower model
#'       value.
#'   (2) CG iterate exits the trust region: truncate to the boundary.
#'   (3) Model gradient is small enough: accept the interior CG iterate.
#'
#' Reference: Nocedal & Wright, "Numerical Optimization" (2006), Algorithm 7.2.
#'
#' @param gn           Gradient vector (length n).
#' @param H            Hessian matrix (n x n); need not be positive definite.
#' @param trust_radius Trust-region radius (positive scalar).
#' @param tol          Convergence tolerance for the model gradient. Defaults to
#'                     min(0.5, sqrt(||g||)) * ||g||, the forcing function from
#'                     Nocedal & Wright that gives superlinear convergence.
#' @param max_iter     Maximum CG iterations (default: length(gn)).
#'
#' @return Numeric vector p of length n, the proposed step, with ||p|| <=
#'         trust_radius.
quadratic_trust_region_update <- function(gn, H, trust_radius,
                                tol = NULL, max_iter = NULL) {
  n      <- length(gn)
  g_norm <- sqrt(sum(gn^2))

  # Default tolerance: "superlinear" forcing function (Nocedal & Wright p. 171)
  if (is.null(tol)){tol      <- min(0.5, sqrt(g_norm)) * g_norm}
  if (is.null(max_iter)){ max_iter <- n }

  # ---------------------------------------------------------------------------
  # Helper: boundary intersection
  #   Find tau satisfying  ||z + tau*d||^2 = trust_radius^2.
  #   Since z is strictly inside the trust region (||z|| < trust_radius),
  #   the quadratic always has one positive and one negative real root.
  #   Returns the positive root by default; both roots when both_roots = TRUE.
  # ---------------------------------------------------------------------------
  boundary_tau <- function(z, d, both_roots = FALSE) {
    a    <- sum(d * d)
    b    <- 2.0 * sum(z * d)
    cc   <- sum(z * z) - trust_radius^2   # < 0 when z is inside the region
    disc <- sqrt(max(b^2 - 4.0 * a * cc, 0.0))
    tau_pos <- (-b + disc) / (2.0 * a)
    if (both_roots){return(c((-b - disc) / (2.0 * a), tau_pos))}
    return(tau_pos)
  }

  # Model value m(p) = g'p + (1/2) p'Hp  (only evaluated for negative-
  # curvature boundary choice, so this path is rarely taken).
  model_val <- function(p) {
    return(sum(gn * p) + 0.5 * as.numeric(crossprod(p, H %*% p)))
  }

  # ---------------------------------------------------------------------------
  # Initialise
  # ---------------------------------------------------------------------------
  z   <- rep(0.0, n)    # current iterate, starts at the origin
  r   <- gn             # r_j = grad m(z_j) = g + H z_j; initially g + H*0 = g
  d   <- -r             # first CG direction = negative gradient
  rtr <- sum(r * r)     # r'r, cached to avoid recomputation

  # Trivial early exit: gradient is already negligibly small
  if (g_norm < tol){return(list(params = z, predicted_red = 0))}

  # ---------------------------------------------------------------------------
  # CG iterations
  # ---------------------------------------------------------------------------
  for (j in seq_len(max_iter)) {

    Hd  <- as.numeric(H %*% d)     # one matrix-vector product per iteration
    kap <- sum(d * Hd)             # curvature  kappa = d'Hd

    # ---- (1) Negative or zero curvature ------------------------------------
    # The quadratic model is concave along d; step to the boundary and pick
    # the intersection point (tau_neg or tau_pos) with the lower model value.
    if (kap <= 0) {
      taus <- boundary_tau(z, d, both_roots = TRUE)
      p1   <- z + taus[1] * d
      p2   <- z + taus[2] * d
      if(model_val(p1) <= model_val(p2)){
        return(list(params = p1, predicted_red = -model_val(p1)))
      } else {
        return(list(params = p2, predicted_red = -model_val(p2)))
      }
    }

    alpha <- rtr / kap       # CG step length
    z_new <- z + alpha * d   # candidate new iterate

    # ---- (2) Step exits the trust region -----------------------------------
    # Truncate to the boundary along the current direction.
    if (sqrt(sum(z_new * z_new)) >= trust_radius) {
      tau <- boundary_tau(z, d)   # positive root only
      return(list(params = (z + tau * d), predicted_red = -model_val(z + tau * d)))
    }

    # Update residual  r_{j+1} = r_j + alpha_j * H d_j = g + H z_{j+1}
    r_new   <- r + alpha * Hd
    rtr_new <- sum(r_new * r_new)

    # ---- (3) Residual small enough: accept interior step -------------------
    if (sqrt(rtr_new) < tol) return(list(params = z_new, predicted_red = -model_val(z_new)))

    # Continue CG
    beta <- rtr_new / rtr
    d    <- -r_new + beta * d
    z    <- z_new
    r    <- r_new
    rtr  <- rtr_new
  }

  # Maximum iterations reached: return the last accepted interior iterate
  return(list(params = z, predicted_red = model_val(z)))
}

#' Update the secant matrix with the SR1 update rule
#' @description
#' See 6.2 in Noccedal and Wright
#' @param Bk N&W
#' @param sk N&W
#' @param r N&W
#' @param yk N&W
#'
#' @export
SR1_update <- function(Bk, sk, r, yk){
  secDiff <- yk - Bk %*% sk
  if(abs(sk%*% secDiff >= r*sum(sk^2)*sum(secDiff^2))){
    return(Bk + (secDiff %*% t(secDiff))/(secDiff %*%  sk))
  } else {
    return(Bk)
  }
}





















