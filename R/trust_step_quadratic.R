#' One trust region step by quadratic interpolation
#' @description
#' An alternative to \code{\link{trust_step_conjugate}}'s sweep of line searches, in the
#' spirit of UOBYQA (Powell 2002): fit a quadratic model to the contrast function on a set
#' of interpolation points, minimise that model over the trust region, and take the
#' minimiser as the step.
#'
#' Interpolation needs the feasible region mapped out before the points are chosen, which
#' is the substantive difference from the line search. A line search can meet a point the
#' importance sampling cannot support, drop it, and carry on --
#' \code{\link{trust_line_evaluate}} leaves `NA` there and `which.min` passes over it. An
#' interpolation set cannot: a missing value is a missing row of the interpolation system,
#' which then has no unique solution. So the step opens by bisecting along each coordinate
#' for the point where `ess/nSims` crosses `eta_trust` (see \code{\link{ess_box_bisect}}),
#' and places its points inside the resulting box.
#'
#' Within one step every evaluation re-weights the same ensemble, so the contrast is a
#' deterministic, smooth function of the parameters and fitting an exact interpolant to it
#' is well posed. The sampling error that \code{\link{trust_region_loop}} warns about lives
#' between ensembles, not within one, and so does not argue for a regression model here.
#'
#' @section The box and its corners:
#' The box is the product of the per-coordinate crossings, and is not contained in the
#' region where the requirement actually holds: loss of effective sample size compounds
#' across coordinates, so a point moving several of them at once is worse supported than
#' any of the single-coordinate moves that make it up. Measured on study1's starting
#' parameters at `nSims = 300`, `delta = 0.2` and `eta_trust = 0.6`, the crossings sit at
#' `(0.086, 0.200, 0.169)` in `(log kappa, log omega, log mu)`, each at `ess/nSims = 0.6`
#' by construction, while the corner those three span comes to 0.409 -- well inside the
#' region the requirement rejects. It takes a uniform shrink to about 0.6 of the corner
#' before it passes.
#'
#' Every interpolation point is therefore checked and, if it fails, pulled back along its
#' own ray from `x_0` to the last position that passes (see \code{\link{ess_ray_limit}}).
#' That keeps the full width of the box on the coordinates that can take it -- which
#' matters, since the region is markedly anisotropic, omega reaching the trust radius
#' without ever crossing while kappa crosses at less than half of it -- without ever
#' leaving a hole in the interpolation set.
#'
#' @section Cost:
#' A step costs `2n` ray searches for the box, at most `bisection_iterations + 1`
#' evaluations each, plus one evaluation per interpolation point and its repairs, plus one
#' at the model minimiser. For three free parameters at the defaults that is around 60 to
#' 70 evaluations against the several hundred a `max.directions = NULL` conjugate step
#' spends, which is the reason to reach for this method.
#'
#' The evaluations are not interchangeable, though, and the difference is large enough to
#' plan around. `trust_function` recomputes the `O(n_daughter x n_parent)` daughter kernel
#' sums whenever omega moves off the value its cache holds, and reuses them when it does
#' not. At `nSims = 300` a full evaluation measured 3.2 ms with omega held and 23.4 ms with
#' omega moved, a factor of 7. Two consequences are built into this method: a ray search
#' tries the trust radius before it bisects, so a coordinate on which `delta` binds before
#' the requirement does costs one evaluation instead of `bisection_iterations + 1` (on the
#' measurement above this is exactly what omega does), and the interpolation points are
#' evaluated grouped by omega, so a set holding three distinct omegas pays for three sets
#' of kernel sums rather than one per point.
#'
#' @param x_0 Parameter values the ensemble behind `trust_function` was simulated from,
#' and the centre of both the trust region and the interpolation set.
#' @param delta Trust region radius.
#' @param trust_function Trust function. Called on a parameter vector and expected to
#' return a list with `f_est` (the contrast estimate) and `ess` (the effective sample size
#' of the importance sampling weights behind it).
#' @param eta_trust Smallest `ess/nSims` an iterate is allowed to have, between 0 and 1.
#' @param nSims Number of simulations the ensemble behind `trust_function` holds.
#' @param interp_fraction How much of the box to spread the interpolation points over, as
#' a fraction of the distance from `x_0` to each of its faces. The model is fitted on that
#' smaller box and then minimised over the whole region, which is the separation UOBYQA
#' keeps between the radius a model is built on and the radius it is trusted over.
#'
#' It matters, and not subtly: a quadratic is a good model of the contrast only close to
#' `x_0`. Fitting one across the full box and minimising it produced a model predicting a
#' *negative* contrast where the true value was 0.19, while the same set at a fraction of
#' 0.15 predicted 0.1295 against a true 0.1295. Over 8 ensembles at `nSims = 300` the
#' median contrast reached was
#'
#' \tabular{lrrrr}{
#'   \tab `0.1` \tab `0.25` \tab `0.5` \tab `1.0` \cr
#'   `delta = 0.05` \tab 0.1329 \tab 0.1330 \tab 0.1288 \tab 0.1296 \cr
#'   `delta = 0.1`  \tab 0.1309 \tab 0.1278 \tab 0.1364 \tab 0.1523 \cr
#'   `delta = 0.2`  \tab 0.1408 \tab 0.1368 \tab 0.1438 \tab 0.1471
#' }
#'
#' which is why the default is 0.25. The best value drifts down as `delta` grows, since a
#' fixed fraction of a larger box is a larger region for the quadratic to have to describe,
#' so a fit run at a large trust radius may do better at a smaller fraction than this.
#' @param bisection_iterations How many bisection steps each ray search takes after it has
#' established that the requirement fails at the far end. The interval starts at the trust
#' radius and halves each time, so the crossing is located to `delta/2^bisection_iterations`
#' -- at the default of 10 and a radius of 0.2, to about `2e-4`, which is far finer than
#' the box needs to be. It is worth lowering before anything else if a step is too slow.
#' @param ess_bounds Per-parameter breakdown bounds to screen against and add to, from
#' \code{\link{ess_bounds_init}}, or `NULL` to start from scratch. Bounds carried over
#' from a previous step at the same `x_0` matter more here than to the conjugate method,
#' since they let the ray searches narrow a bracket without paying for the probes above a
#' crossing this centre has already located. Bounds belonging to any other centre are
#' discarded. See \code{\link{trust_region_loop}}'s `carry.ess_bounds`.
#'
#' @return A list with `x_star` (the step taken), `f_pred` and `ess_star` (the contrast
#' estimate and effective sample size there), `ess_bounds`, `directions` (always `NULL`,
#' this method keeping no direction set; the element is present so that the two methods
#' return the same shape), and `box` (the per-coordinate `lower` and `upper` offsets from
#' `x_0` that the interpolation was carried out in).
#'
#' @export
trust_step_quadratic <- function(x_0, delta, trust_function, eta_trust, nSims,
                                 interp_fraction = 0.25, bisection_iterations = 10,
                                 ess_bounds = NULL){
  if((interp_fraction <= 0) || (interp_fraction > 1)){
    stop("interp_fraction must be in (0, 1].")
  }
  n <- length(x_0)
  state <- list(ess_bounds = ess_bounds_carried(ess_bounds, x_0), evals = list())

  # 1. Map the region ----
  boxRes <- ess_box_bisect(x_0 = x_0, delta = delta, trust_function = trust_function,
                           eta_trust = eta_trust, nSims = nSims,
                           bisection_iterations = bisection_iterations, state = state)
  state <- boxRes$state
  box <- boxRes$box

  # 2. Place and evaluate the interpolation set ----
  # Fitted on a smaller box than it is minimised over: see `interp_fraction`.
  interp_box <- list(lower = box$lower*interp_fraction, upper = box$upper*interp_fraction)
  setRes <- quadratic_interp_set(x_0 = x_0, box = interp_box, delta = delta,
                                 trust_function = trust_function, eta_trust = eta_trust,
                                 nSims = nSims, bisection_iterations = bisection_iterations,
                                 state = state)
  state <- setRes$state

  # 3. Fit the model and minimise it ----
  # A set that could not be completed leaves the interpolation system without a unique
  # solution. Rather than fitting something arbitrary to it, the step falls back on the
  # best point the mapping and the set have already paid for, which is a genuine
  # improvement whenever one was found and x_0 itself otherwise.
  x_star <- trust_best_eval(state$evals, default = x_0)
  model <- quadratic_model_fit(S = setRes$S, f = setRes$f, n = n)
  if(!is.null(model)){
    s_model <- quadratic_model_minimise(model = model, lower = box$lower,
                                        upper = box$upper, delta = delta,
                                        starts = quadratic_model_starts(setRes$S, box))
    if(any(s_model != 0)){
      # The minimiser is a prediction, and the box it was found in over-states the region,
      # so it is walked back along its own ray until it is a point the ensemble supports.
      rayRes <- ess_ray_limit(x_0 = x_0, dir = s_model, t_hi = 1,
                              trust_function = trust_function, eta_trust = eta_trust,
                              nSims = nSims, bisection_iterations = bisection_iterations,
                              state = state)
      state <- rayRes$state
      # Whatever the model predicted, the step returned is the best point actually
      # evaluated. The model chooses where to look; it does not get to overrule a
      # measurement taken at the point it chose.
      x_star <- trust_best_eval(state$evals, default = x_star)
    }
  }

  starRes <- trust_function(x_star)
  return(list(x_star = x_star,
              f_pred = starRes$f_est,
              ess_star = starRes$ess,
              ess_bounds = state$ess_bounds,
              directions = NULL,
              box = box))
}

#' Evaluate the trust function at one point and record the outcome
#' @description
#' The single place this method calls `trust_function`. A point the recorded bounds
#' already rule out is not evaluated at all; a point that is evaluated has its outcome
#' written into the bounds either way, and is kept in the running list of evaluations so
#' that the step can fall back on the best point it has actually measured.
#'
#' @param x Parameter vector to evaluate.
#' @param trust_function Trust function, returning a list with `f_est` and `ess`.
#' @param eta_trust Smallest `ess/nSims` an iterate is allowed to have.
#' @param nSims Number of simulations the ensemble behind `trust_function` holds.
#' @param state A list with `ess_bounds` and `evals`, threaded through the step.
#'
#' @return A list with `passed`, `f` (`NA` unless the point passed), and the updated
#' `state`.
trust_eval_record <- function(x, trust_function, eta_trust, nSims, state){
  if(!ess_bounds_feasible(state$ess_bounds, x)){
    return(list(passed = FALSE, f = NA_real_, state = state))
  }
  res <- trust_function(x)
  passed <- (res$ess/nSims > eta_trust) && is.finite(res$f_est)
  state$ess_bounds <- ess_bounds_record(state$ess_bounds, x, passed = passed)
  if(passed){
    state$evals[[length(state$evals) + 1]] <- list(x = x, f = res$f_est)
  }
  return(list(passed = passed, f = if(passed) res$f_est else NA_real_, state = state))
}

#' Best point evaluated so far
#' @description
#' Picks the smallest contrast value out of the evaluations recorded by
#' \code{\link{trust_eval_record}}. Only points that passed the effective sample size
#' requirement are ever recorded, so anything in the list is a legitimate step.
#'
#' @param evals The `evals` element of the step's state.
#' @param default Point to return when nothing has been recorded.
#'
#' @return A parameter vector.
trust_best_eval <- function(evals, default){
  if(length(evals) == 0){
    return(default)
  }
  f <- vapply(evals, function(e) e$f, numeric(1))
  return(evals[[which.min(f)]]$x)
}

#' Largest feasible step along a ray
#' @description
#' Finds how far `x_0 + t*dir` can go, for `t` in `(0, t_hi]`, before `ess/nSims` drops to
#' `eta_trust`. The far end is tried first, so a ray that is feasible all the way to
#' `t_hi` costs one evaluation rather than a full bisection -- which is the common case on
#' a coordinate where the trust radius binds before the importance sampling does.
#'
#' Bisection is the right search here because the effective sample size falls off
#' monotonically as a point moves away from the parameters the ensemble was simulated
#' from, so the ray crosses the requirement exactly once. Measured on study1's starting
#' parameters at `nSims = 300`, `ess/nSims` on a 13-point grid across `+/- 0.2` was
#' monotone decreasing outwards on both sides of all three coordinate axes, with a single
#' crossing per side.
#'
#' @param x_0 Centre of the trust region, which always passes: every importance sampling
#' weight against its own ensemble is 1, so the effective sample size there is `nSims`.
#' @param dir Direction of the ray. Not required to be a unit vector; `t` is in units of
#' `dir`.
#' @param t_hi Largest `t` to consider.
#' @param trust_function Trust function, returning a list with `f_est` and `ess`.
#' @param eta_trust Smallest `ess/nSims` an iterate is allowed to have.
#' @param nSims Number of simulations the ensemble behind `trust_function` holds.
#' @param bisection_iterations Bisection steps to take once the far end has failed.
#' @param state A list with `ess_bounds` and `evals`, threaded through the step.
#'
#' @return A list with `t` (the largest value found to pass, 0 if none did) and the
#' updated `state`.
ess_ray_limit <- function(x_0, dir, t_hi, trust_function, eta_trust, nSims,
                          bisection_iterations, state){
  hiRes <- trust_eval_record(x = x_0 + t_hi*dir, trust_function = trust_function,
                             eta_trust = eta_trust, nSims = nSims, state = state)
  state <- hiRes$state
  if(hiRes$passed){
    return(list(t = t_hi, state = state))
  }

  # x_0 passes and t_hi does not, so the crossing is bracketed. `lo` is only ever moved to
  # a value that has been evaluated and passed, so what is returned is a point the
  # ensemble supports rather than an estimate of where support ends.
  lo <- 0
  hi <- t_hi
  for(k in seq_len(bisection_iterations)){
    mid <- (lo + hi)/2
    midRes <- trust_eval_record(x = x_0 + mid*dir, trust_function = trust_function,
                                eta_trust = eta_trust, nSims = nSims, state = state)
    state <- midRes$state
    if(midRes$passed){
      lo <- mid
    } else {
      hi <- mid
    }
  }
  return(list(t = lo, state = state))
}

#' Per-coordinate box the importance sampling supports
#' @description
#' Runs \code{\link{ess_ray_limit}} out along each coordinate axis in both directions,
#' capped at the trust radius, and returns the offsets from `x_0` it reached. The box is
#' where \code{\link{trust_step_quadratic}} places its interpolation points; see that
#' function on why its corners still have to be checked individually.
#'
#' @param x_0 Centre of the trust region.
#' @param delta Trust region radius, which caps each axis: along a single coordinate the
#' trust region reaches exactly `delta`.
#' @param trust_function Trust function, returning a list with `f_est` and `ess`.
#' @param eta_trust Smallest `ess/nSims` an iterate is allowed to have.
#' @param nSims Number of simulations the ensemble behind `trust_function` holds.
#' @param bisection_iterations Bisection steps per ray search.
#' @param state A list with `ess_bounds` and `evals`, threaded through the step.
#'
#' @return A list with `box` (a list of per-coordinate `lower` and `upper` offsets from
#' `x_0`, `lower <= 0 <= upper`) and the updated `state`.
ess_box_bisect <- function(x_0, delta, trust_function, eta_trust, nSims,
                           bisection_iterations, state){
  n <- length(x_0)
  lower <- rep(0, n)
  upper <- rep(0, n)
  for(i in seq_len(n)){
    e_i <- rep(0, n)
    e_i[i] <- 1
    upRes <- ess_ray_limit(x_0 = x_0, dir = e_i, t_hi = delta,
                           trust_function = trust_function, eta_trust = eta_trust,
                           nSims = nSims, bisection_iterations = bisection_iterations,
                           state = state)
    state <- upRes$state
    upper[i] <- upRes$t
    downRes <- ess_ray_limit(x_0 = x_0, dir = -e_i, t_hi = delta,
                             trust_function = trust_function, eta_trust = eta_trust,
                             nSims = nSims, bisection_iterations = bisection_iterations,
                             state = state)
    state <- downRes$state
    lower[i] <- -downRes$t
  }
  return(list(box = list(lower = lower, upper = upper), state = state))
}

#' Interpolation points for a quadratic model, and the contrast at each
#' @description
#' Places the `(n+1)(n+2)/2` points a quadratic in `n` variables is determined by: `x_0`
#' itself, one point out along each coordinate in each direction, and one point per pair
#' of coordinates moving both at once. The axial points fix the constant, the gradient and
#' the diagonal of the Hessian; the pair points fix its off-diagonal entries.
#'
#' Points are taken to the edge of the box, so the model is fitted over the whole region it
#' will then be minimised in. Any point the ensemble turns out not to support -- in
#' practice the pair points, which is exactly where the box over-states the region -- is
#' pulled back along its own ray from `x_0` until it passes, rather than dropped, since a
#' dropped point leaves the interpolation system underdetermined.
#'
#' The points are visited in order of their omega coordinate. `trust_function` keeps the
#' daughter kernel sums for one omega at a time and recomputing them costs about seven
#' times an evaluation that does not, so the points that share an omega -- the two axial
#' points of every other coordinate, in particular -- are cheaper taken together than
#' interleaved with points that move it. The saving is partial rather than total, since a
#' point that needs repair walks its own ray and moves omega again on the way.
#'
#' @param x_0 Centre of the trust region and of the interpolation set.
#' @param box Per-coordinate offsets from \code{\link{ess_box_bisect}}.
#' @param delta Trust region radius.
#' @param trust_function Trust function, returning a list with `f_est` and `ess`.
#' @param eta_trust Smallest `ess/nSims` an iterate is allowed to have.
#' @param nSims Number of simulations the ensemble behind `trust_function` holds.
#' @param bisection_iterations Bisection steps per repair.
#' @param state A list with `ess_bounds` and `evals`, threaded through the step.
#'
#' @return A list with `S` (the offsets from `x_0` actually used, one row per point),
#' `f` (the contrast at each, `NA` where no feasible position was found) and the updated
#' `state`.
quadratic_interp_set <- function(x_0, box, delta, trust_function, eta_trust, nSims,
                                 bisection_iterations, state){
  n <- length(x_0)

  # Offsets first, evaluations second, so the points can be reordered by omega before any
  # of them is paid for.
  S <- matrix(0, nrow = 1, ncol = n)
  for(i in seq_len(n)){
    for(side in c(box$upper[i], box$lower[i])){
      s <- rep(0, n)
      s[i] <- side
      S <- rbind(S, s)
    }
  }
  if(n > 1){
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        s <- rep(0, n)
        # The wider side of each coordinate, so the pair point carries as much
        # information about the cross term as the region allows before repair.
        s[i] <- if(box$upper[i] >= -box$lower[i]) box$upper[i] else box$lower[i]
        s[j] <- if(box$upper[j] >= -box$lower[j]) box$upper[j] else box$lower[j]
        S <- rbind(S, s)
      }
    }
  }
  rownames(S) <- NULL

  f <- rep(NA_real_, nrow(S))
  # Free coordinate 2 is omega whenever the free set starts at log kappa, which covers
  # both the main fit and the K-function block. On any other subset this groups by
  # whatever coordinate happens to sit there, which costs nothing beyond the ordering
  # being no better than arbitrary.
  omega_col <- if(n >= 2) 2 else 1
  visit <- order(S[, omega_col])
  for(idx in visit){
    if(all(S[idx, ] == 0)){
      # x_0 itself. Always feasible, and evaluating it is how the model gets its constant.
      res <- trust_eval_record(x = x_0, trust_function = trust_function,
                               eta_trust = eta_trust, nSims = nSims, state = state)
      state <- res$state
      f[idx] <- res$f
      next
    }
    rayRes <- ess_ray_limit(x_0 = x_0, dir = S[idx, ], t_hi = 1,
                            trust_function = trust_function, eta_trust = eta_trust,
                            nSims = nSims, bisection_iterations = bisection_iterations,
                            state = state)
    state <- rayRes$state
    S[idx, ] <- rayRes$t * S[idx, ]
    if(rayRes$t > 0){
      f[idx] <- trust_eval_lookup(state$evals, x_0 + S[idx, ])
    }
  }

  return(list(S = S, f = f, state = state))
}

#' Contrast value already recorded at a point
#' @description
#' \code{\link{ess_ray_limit}} evaluates the position it returns on the way to finding it,
#' so the value is already in the step's evaluation list and does not need recomputing.
#'
#' @param evals The `evals` element of the step's state.
#' @param x Parameter vector to look up.
#'
#' @return The contrast value, or `NA` if the point is not in the list.
trust_eval_lookup <- function(evals, x){
  for(e in rev(evals)){
    if(isTRUE(all.equal(e$x, x, tolerance = 0))){
      return(e$f)
    }
  }
  return(NA_real_)
}

#' Fit a quadratic to the interpolation set
#' @description
#' Solves for the `c`, `g` and `H` of `m(s) = c + g's + 0.5 s'Hs` that reproduce the
#' contrast at every interpolation point. The basis is ordered constant, then the linear
#' terms, then the diagonal quadratic terms, then the off-diagonal ones, which is `1 + 2n
#' + n(n-1)/2 = (n+1)(n+2)/2` coefficients for as many points.
#'
#' The system is solved in the least squares sense, so that a set which repair has left
#' slightly over- or under-determined still produces a model rather than an error. A set
#' that has lost a dimension outright -- repair having collapsed points onto each other or
#' onto `x_0` -- is rank deficient, and that is reported as `NULL` rather than fitted,
#' since a model fitted through a degenerate set would send the step somewhere its points
#' say nothing about.
#'
#' @param S Interpolation offsets from `x_0`, one row per point.
#' @param f Contrast value at each point, possibly with `NA`s.
#' @param n Number of free parameters.
#'
#' @return A list with `g` and `H`, or `NULL` when the set will not support a model.
quadratic_model_fit <- function(S, f, n){
  keep <- is.finite(f)
  if(sum(keep) < (n + 1)*(n + 2)/2){
    return(NULL)
  }
  S <- S[keep, , drop = FALSE]
  f <- f[keep]

  basis <- function(s){
    cross <- numeric(0)
    if(n > 1){
      for(i in 1:(n-1)){
        for(j in (i+1):n){
          cross <- c(cross, s[i]*s[j])
        }
      }
    }
    return(c(1, s, 0.5*s^2, cross))
  }
  A <- t(apply(S, 1, basis))
  qrA <- qr(A)
  if(qrA$rank < ncol(A)){
    return(NULL)
  }
  coef <- qr.coef(qrA, f)
  if(any(!is.finite(coef))){
    return(NULL)
  }

  g <- coef[2:(n + 1)]
  H <- matrix(0, nrow = n, ncol = n)
  diag(H) <- coef[(n + 2):(2*n + 1)]
  if(n > 1){
    k <- 2*n + 2
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        H[i, j] <- coef[k]
        H[j, i] <- coef[k]
        k <- k + 1
      }
    }
  }
  return(list(g = as.numeric(g), H = H))
}

#' Starting points for the model minimisation
#' @description
#' The model is a general quadratic and its Hessian is whatever the contrast happened to
#' look like over the interpolation set, so it may well be indefinite and have several
#' local minima on the region. Projected gradient descent is therefore started from the
#' centre, from each vertex of the box, and from each interpolation point, and the best
#' result is kept.
#'
#' @param S Interpolation offsets actually used.
#' @param box Per-coordinate offsets from \code{\link{ess_box_bisect}}.
#'
#' @return A matrix of starting offsets, one row each.
quadratic_model_starts <- function(S, box){
  n <- length(box$lower)
  starts <- rbind(rep(0, n), S)
  for(i in seq_len(n)){
    s <- rep(0, n)
    s[i] <- box$upper[i]
    starts <- rbind(starts, s)
    s[i] <- box$lower[i]
    starts <- rbind(starts, s)
  }
  return(unique(starts))
}

#' Minimise the quadratic model over the trust region
#' @description
#' Projected gradient descent with backtracking, over the intersection of the box and the
#' ball of radius `delta`. Both sets are convex and both contain the centre, and a point is
#' returned to the intersection by scaling it onto the ball and then clipping it to the
#' box -- clipping only ever moves a coordinate towards zero, so it cannot push the point
#' back outside the ball, and one pass of the two therefore lands inside both.
#'
#' @param model The `g` and `H` from \code{\link{quadratic_model_fit}}.
#' @param lower,upper Per-coordinate offsets bounding the box.
#' @param delta Trust region radius.
#' @param starts Starting offsets, one row each.
#' @param iterations Projected gradient steps per start.
#'
#' @return The offset from `x_0` minimising the model over the region.
quadratic_model_minimise <- function(model, lower, upper, delta, starts, iterations = 200){
  g <- model$g
  H <- model$H
  m <- function(s) sum(g*s) + 0.5*sum(s*as.vector(H %*% s))
  grad <- function(s) g + as.vector(H %*% s)
  project <- function(s){
    nrm <- sqrt(sum(s^2))
    if(nrm > delta){
      s <- s*(delta/nrm)
    }
    return(pmin(pmax(s, lower), upper))
  }

  best_s <- rep(0, length(g))
  best_m <- 0
  for(r in seq_len(nrow(starts))){
    s <- project(starts[r, ])
    val <- m(s)
    for(k in seq_len(iterations)){
      d <- grad(s)
      dNorm <- sqrt(sum(d^2))
      if(!is.finite(dNorm) || dNorm == 0){
        break
      }
      # Backtracking on the projected step. The projection can shorten the step, so the
      # decrease is tested on the projected point rather than the unprojected one.
      stepped <- FALSE
      t <- delta/dNorm
      for(b in 1:30){
        s_new <- project(s - t*d)
        val_new <- m(s_new)
        if(is.finite(val_new) && (val_new < val - 1e-14)){
          s <- s_new
          val <- val_new
          stepped <- TRUE
          break
        }
        t <- t/2
      }
      if(!stepped){
        break
      }
    }
    if(val < best_m){
      best_m <- val
      best_s <- s
    }
  }
  return(best_s)
}
