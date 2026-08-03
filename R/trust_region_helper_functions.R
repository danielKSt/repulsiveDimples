

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
trust_region_update <- function(gn, H, trust_radius,
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
