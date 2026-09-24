#' Do one iteration within trust region
#' @description
#' Minimises the importance sampling estimate of the contrast function over the trust
#' region, by one of two methods, and returns the minimiser as the step to take. This is
#' the dispatcher; \code{\link{trust_step_conjugate}} and
#' \code{\link{trust_step_quadratic}} are the methods themselves and document what each
#' one does.
#'
#' Both re-weight the single ensemble simulated at `x_0`, so both are restricted to the
#' ball of radius `delta` around it and to points where the importance sampling still has
#' an effective sample size of at least `eta_trust*nSims`. What differs is how they search
#' and, more to the point, how they treat a point that fails that requirement.
#'
#' `"conjugate"` sweeps a set of line searches in the manner of Powell (1964), building
#' conjugate directions as it goes. A point it cannot support is simply dropped from the
#' line it was found on and the search carries on, so the requirement costs it nothing
#' beyond the evaluation. This is the original method, and was the default until
#' 2026-09-24.
#'
#' `"quadratic"` fits a quadratic model to the contrast on a set of interpolation points
#' and minimises that, in the manner of UOBYQA (Powell 2002). It cannot drop a point --
#' a missing value is a missing row of the interpolation system -- so it first maps the
#' region the ensemble supports by bisecting along each coordinate, and places its points
#' inside what it finds. It spends far fewer evaluations than a conjugate step with
#' `max.directions = NULL`, and unlike a line search it uses the curvature of the contrast
#' rather than only its values along a line. This is the default.
#'
#' Measured on study1 with 200 paired fits at each of four ensemble sizes -- same seeds,
#' same starting point, same pre-fitting -- the quadratic step is 4.1 to 5.1 times
#' cheaper in CPU across the range, and 5.0 times over the whole study. Its accuracy
#' matches or beats the conjugate step up to `nSims = 1000`: at 500 it is better and much
#' less erratic (spread 0.138 against 0.242, where one conjugate fit ended more than 1.0
#' out in log space). Above that it falls behind, significantly at 10000, where the median
#' error in log parameters is 0.044 against 0.030. That gap is variance rather than bias:
#' at 5000 and 10000 the quadratic step is actually the less biased of the two (0.0125
#' against 0.0243 at 10000), but its spread stops shrinking between those sizes while the
#' conjugate step's keeps falling.
#'
#' Compared at equal cost rather than equal `nSims`, which is the fair comparison when one
#' method is five times cheaper, the quadratic step is ahead for any budget up to about
#' 2000 s a fit -- the conjugate step at 1000 costs what the quadratic step does at 5000,
#' and reaches 0.074 against 0.046. That is why it is the default. The conjugate step is
#' worth its price only when the last of the accuracy matters and five times the compute
#' is affordable; because the quadratic step's spread plateaus, it cannot reach the
#' conjugate step's best at any ensemble size.
#'
#' The two also fail differently. The conjugate step degrades gracefully into a cruder
#' search when the region is awkward, while the quadratic step depends on an
#' interpolation set that the region has not deformed too badly, and falls back on the
#' best point it evaluated when it has.
#'
#' @param x_0 Parameter values used to simulate patternSim
#' @param delta Trust region radius
#' @param trust_function Trust function. Called on a parameter vector and expected to
#' return a list with `f_est` (the contrast estimate) and `ess` (the effective sample
#' size of the importance sampling weights behind it).
#' @param eta_trust Smallest `ess/nSims` an iterate is allowed to have, between 0 and 1.
#' @param nSims Number of simulations the ensemble behind `trust_function` holds.
#' @param method Which method to use, `"quadratic"` (the default) or `"conjugate"`.
#' @param subsection_count,line_iterations,max.directions,directions Passed to
#' \code{\link{trust_step_conjugate}}, and ignored by the quadratic method.
#' @param interp_fraction,bisection_iterations Passed to
#' \code{\link{trust_step_quadratic}}, and ignored by the conjugate method.
#' @param ess_bounds Breakdown bounds from a previous step at the same `x_0`, used by both
#' methods to skip points already known to break the effective sample size requirement, or
#' `NULL` to start from scratch.
#'
#' @return A list with `x_star` (the trust region minimiser), `f_pred` and `ess_star` (the
#' contrast estimate and effective sample size there), `ess_bounds` (the per-parameter
#' breakdown bounds accumulated over this step), and `directions` (the direction set as
#' this step left it, or `NULL` from the quadratic method, which keeps none). The
#' quadratic method additionally returns `box`.
#'
#' @export
trust_step <- function(x_0, delta, trust_function, eta_trust, nSims,
                       method = c("quadratic", "conjugate"),
                       subsection_count = 11, line_iterations = 4, max.directions = 1,
                       directions = NULL, interp_fraction = 0.25,
                       bisection_iterations = 10, ess_bounds = NULL){
  method <- match.arg(method)
  if(method == "quadratic"){
    return(trust_step_quadratic(x_0 = x_0, delta = delta, trust_function = trust_function,
                                eta_trust = eta_trust, nSims = nSims,
                                interp_fraction = interp_fraction,
                                bisection_iterations = bisection_iterations,
                                ess_bounds = ess_bounds))
  }
  return(trust_step_conjugate(x_0 = x_0, delta = delta, trust_function = trust_function,
                              eta_trust = eta_trust, nSims = nSims,
                              subsection_count = subsection_count,
                              line_iterations = line_iterations,
                              max.directions = max.directions,
                              directions = directions,
                              ess_bounds = ess_bounds))
}

#' One trust region step by conjugate direction line searches
#' @description
#' Iterate from params_sim using importance sampling on the realizations patternSim, simulated from params_sim.
#' The iteration is restricted to the ball centered on the parameters used for simulation,
#' with radius equal to the trust region radius.
#'
#' Every evaluation re-weights the one ensemble simulated at `x_0`, so it is only worth
#' trusting while the importance sampling has not degenerated. Iterates are therefore
#' additionally required to satisfy `ess/nSims > eta_trust`, and the parameter values at
#' which that requirement was seen to break down are carried between the line searches
#' (see \code{\link{ess_bounds_init}}) so that later ones can skip points already known
#' to be too far out.
#'
#' @param x_0 Parameter values used to simulate patternSim
#' @param delta Trust region radius
#' @param trust_function Trust function. Called on a parameter vector and expected to
#' return a list with `f_est` (the contrast estimate) and `ess` (the effective sample
#' size of the importance sampling weights behind it).
#' @param eta_trust Smallest `ess/nSims` an iterate is allowed to have, between 0 and 1.
#' @param nSims Number of simulations the ensemble behind `trust_function` holds, i.e.
#' the largest effective sample size attainable.
#' @param subsection_count Points per line search grid, at least 5. Each refinement
#' re-grids the four intervals around the current minimum, so the search narrows by
#' `4/(subsection_count - 1)` per round -- 0.4 at the default 11, but only 0.8 at 6, and
#' nothing at all at 5.
#' @param line_iterations How many times each line search refines its grid, passed to
#' \code{\link{trust_line_steps}} as `iterations`.
#' @param directions Direction set to start this step from, as an `n x n` matrix of
#' columns, or `NULL` for the coordinate axes. Pass back the `directions` a previous step
#' returned to continue accumulating conjugate directions across trust region iterations
#' instead of restarting from the axes every time -- the accumulation is the mechanism
#' behind Powell (1964), and restarting discards it. Restarting also makes every step open
#' with a line search along free parameter 1, which biases the search toward that
#' parameter. See \code{\link{trust_region_loop}}'s `carry.directions`.
#' @param ess_bounds Breakdown bounds to screen against and add to, from
#' \code{\link{ess_bounds_init}}, or `NULL` to start from scratch. Bounds carried over
#' from a previous step at the same `x_0` let this one skip points already known to break
#' the requirement; bounds belonging to any other centre are discarded. See
#' \code{\link{trust_region_loop}}'s `carry.ess_bounds`.
#' @param max.directions How many times the conjugate direction set may be updated.
#' `NULL` restores the original `2*length(x_0) - 1`. The default of 1 keeps only the
#' first update: a step is then a sweep along the current directions plus one conjugate
#' correction. That makes each step cruder but several times cheaper, and over a whole fit
#' the outer loop compensates by taking more steps -- on study1 the saving was large and
#' unambiguous while the difference in final contrast stayed inside the run-to-run spread,
#' which is not the same as showing the two are equivalent. At 0 the step degenerates into
#' a single line search along the first free parameter.
#'
#' @section Cost:
#' A step costs at most `(1 + max.directions*(n + 1)) * (line_iterations + 1) *
#' subsection_count` calls to `trust_function` for `n` free parameters, before the ESS
#' screening prunes any of them and before the boundary test ends the loop early. At the
#' defaults with three free parameters that is 5 line searches of 55 points each; passing
#' `max.directions = NULL` restores the former 25. `trust_step` dominates the wall time of
#' a fit -- on study1 one `simulation_step` is 0.16 s against 41 s for a `max.directions =
#' NULL` step -- so these three arguments are the levers on how long a fit takes.
#'
#' @return A list with `x_star` (the trust region minimiser), `f_pred` and `ess_star`
#' (the contrast estimate and effective sample size there), `ess_bounds` (the
#' per-parameter breakdown bounds accumulated over this step), and `directions` (the
#' direction set as this step left it, to feed back in through `directions`).
#'
#' @export
trust_step_conjugate <- function(x_0, delta, trust_function, eta_trust, nSims,
                       subsection_count = 11, line_iterations = 4, max.directions = 1,
                       directions = NULL, ess_bounds = NULL){
  n <- length(x_0)
  if(subsection_count < 5){
    stop("subsection_count must be at least 5, since each refinement brackets the grid minimum two points out on either side.")
  }
  # NULL asks for the original cap. The loop below counts updates from zero, so the
  # `k < 2*n` it used to run starting from k = 1 is 2*n - 1 of them.
  if(is.null(max.directions)){
    max.directions <- 2*n - 1
  }
  ess_bounds <- ess_bounds_carried(ess_bounds, x_0)

  if(n == 1){
    lineRes <- trust_line_steps(p_k = c(1), x_k = x_0, x_0 = x_0, delta = delta,
                                trust_function = trust_function, ess_bounds = ess_bounds,
                                eta_trust = eta_trust, nSims = nSims,
                                # A one-dimensional step is a single line search, so it
                                # can afford a finer grid than one of the many that make
                                # up a step over several parameters. The +2 keeps the
                                # 11/13 pairing the defaults have always had.
                                subsection_count = subsection_count + 2,
                                iterations = line_iterations)
    x_k <- lineRes$x
    ess_bounds <- lineRes$ess_bounds
    p <- matrix(data = 1, nrow = 1, ncol = 1)
  } else {
    p <- direction_set_init(directions, n)
    k <- 0
    lineRes <- trust_line_steps(x_k = x_0, p_k = p[, 1], x_0 = x_0, delta = delta,
                                trust_function = trust_function, ess_bounds = ess_bounds,
                                eta_trust = eta_trust, nSims = nSims,
                                subsection_count = subsection_count,
                                iterations = line_iterations)
    x_k <- lineRes$x
    ess_bounds <- lineRes$ess_bounds
    hit_boundary <- (sqrt(sum(abs(x_k-x_0)^2)) > 0.99*delta)
    rotated <- FALSE
    while(!hit_boundary && k < max.directions){
      z_j <- matrix(data = 0, nrow = n, ncol = n + 1)
      z_j[, 1] <- x_k
      for(j in 1:n){
        lineRes <- trust_line_steps(p_k = p[, j], x_k = z_j[, j], x_0 = x_0,
                                    delta = delta, trust_function = trust_function,
                                    ess_bounds = ess_bounds, eta_trust = eta_trust,
                                    nSims = nSims, subsection_count = subsection_count,
                                    iterations = line_iterations)
        z_j[, j+1] <- lineRes$x
        ess_bounds <- lineRes$ess_bounds
      }
      for(j in 1:(n-1)){
        p[, j] <- p[, j+1]
      }
      # Stored as a unit vector. A line search is invariant to the length of its
      # direction -- find_alpha_range rescales alpha to match -- so this moves no
      # evaluated point, but it makes det(p) a scale-free measure of how close the set is
      # to losing a dimension. That matters once the set is carried across iterations
      # rather than rebuilt from the axes each time, since near-degeneracy accumulates.
      v <- z_j[, n + 1] - z_j[, 1]
      vNorm <- sqrt(sum(v^2))
      p[, n] <- if(vNorm > 0) v/vNorm else v
      # Exact test, as before. Note that normalising above changes what this catches: with
      # an unnormalised direction a tiny displacement gives column entries around 1e-10,
      # whose det underflows to exactly zero and fires the rotation fallback even though
      # the directions are not actually dependent. Normalised, only genuine degeneracy
      # trips it. So this differs from the pre-2026-09-14 package in that narrow case, and
      # deliberately -- the old trigger was a floating point artefact. A carried set is
      # additionally checked against a tolerance on the way in, by direction_set_init.
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
      lineRes <- trust_line_steps(x_k = z_j[, n + 1], p_k = p[, n], x_0 = x_0,
                                  delta = delta, trust_function = trust_function,
                                  ess_bounds = ess_bounds, eta_trust = eta_trust,
                                  nSims = nSims, subsection_count = subsection_count,
                                  iterations = line_iterations)
      x_k <- lineRes$x
      ess_bounds <- lineRes$ess_bounds
      k <- k + 1
      hit_boundary <- (sqrt(sum(abs(x_k-x_0)^2)) > 0.99*delta)
    }
  }

  starRes <- trust_function(x_k)
  return(list(x_star = x_k,
              f_pred = starRes$f_est,
              ess_star = starRes$ess,
              ess_bounds = ess_bounds,
              directions = p))
}

#' Direction set to open a trust step with
#' @description
#' Returns `directions` when it is a usable `n x n` set, and the coordinate axes
#' otherwise. A carried set is rejected if it has the wrong shape (the number of free
#' parameters changed), holds anything non-finite, or has gone singular, since a singular
#' set has lost a search dimension and would never recover it on its own.
#'
#' @param directions Direction set from a previous step, or `NULL`.
#' @param n Number of free parameters.
#'
#' @return An `n x n` matrix of direction columns.
direction_set_init <- function(directions, n){
  usable <- !is.null(directions) &&
    is.matrix(directions) && all(dim(directions) == c(n, n)) &&
    all(is.finite(directions)) && (abs(det(directions)) >= 1e-8)
  if(usable){
    return(directions)
  }
  p <- matrix(data = 0, nrow = n, ncol = n)
  diag(p) <- 1
  return(p)
}

#' Line minimization using interpolation with maximum reach
#' @description
#' Line optimization for the trust steps, restricted to points where the importance
#' sampling behind `trust_function` still has an effective sample size of at least
#' `eta_trust*nSims`.
#'
#' Points failing that requirement are recorded in `ess_bounds` and are given no
#' contrast value, so they can neither be returned nor bracketed for refinement, and
#' every later point that the recorded bounds already rule out is skipped without
#' being evaluated at all. `x_k` is assumed to satisfy the requirement itself -- it is
#' either the centre of the trust region or the result of an earlier line search -- and
#' is returned unchanged if the whole line is ruled out.
#'
#' @param p_k Search direction
#' @param x_k Search starting point
#' @param x_0 Center of trust region
#' @param delta Radius of trust region
#' @param subsection_count Number of points on the line search grid, at least 5.
#' @param trust_function Trust function. Called on a parameter vector and expected to
#' return a list with `f_est` and `ess`.
#' @param ess_bounds Per-parameter breakdown bounds to screen against and add to, from
#' \code{\link{ess_bounds_init}}.
#' @param eta_trust Smallest `ess/nSims` an iterate is allowed to have, between 0 and 1.
#' @param nSims Number of simulations the ensemble behind `trust_function` holds.
#' @param iterations How many times to re-grid around the current minimum. Each round
#' brackets it two points out on either side, so the searched span shrinks by
#' `4/(subsection_count - 1)`.
#'
#' @return A list with `x` (the best point found on the line, or `x_k` if none passed the
#' ESS requirement) and `ess_bounds` (the bounds with any breakdown found here added).
#'
#' @export
trust_line_steps <- function(p_k, x_k, x_0, delta, trust_function, ess_bounds,
                             eta_trust, nSims, subsection_count = 11, iterations = 4){
  alpha_range <- find_alpha_range(x_k, x_0, p_k, delta)

  # A zero-width range means the line only meets the trust region at x_k itself, p_k
  # being tangent to the boundary there. There is nothing to search, and without this the
  # grid would be subsection_count copies of x_k, each paid for in full.
  if(alpha_range[1] == alpha_range[2]){
    return(list(x = x_k, ess_bounds = ess_bounds))
  }
  alpha_k <- seq(from = alpha_range[1], to = alpha_range[2], length.out = subsection_count)

  x_star <- x_k
  f_star <- NA_real_
  for(j in 1:(iterations + 1)){
    lineRes <- trust_line_evaluate(alpha_k = alpha_k, p_k = p_k, x_k = x_k,
                                   trust_function = trust_function, ess_bounds = ess_bounds,
                                   eta_trust = eta_trust, nSims = nSims)
    ess_bounds <- lineRes$ess_bounds

    # which.min() passes over the NAs left behind by the ESS requirement, so an empty
    # j_opt means every point on this stretch of the line was ruled out and there is
    # nothing left to bracket. Whatever was found on an earlier, wider grid still stands.
    j_opt <- which.min(lineRes$f_vals)
    if(length(j_opt) == 0){
      break
    }
    if(is.na(f_star) || (lineRes$f_vals[j_opt] < f_star)){
      f_star <- lineRes$f_vals[j_opt]
      x_star <- lineRes$z_k[, j_opt]
    }

    if(j_opt < 3){
      j_opt <- 3
    } else if(j_opt > subsection_count - 2){
      j_opt <- subsection_count - 2
    }
    alpha_k <- seq(from = alpha_k[j_opt-2], to = alpha_k[j_opt+2], length.out = subsection_count)
  }

  return(list(x = x_star, ess_bounds = ess_bounds))
}

#' Evaluate the trust function along one grid of line search steps
#' @description
#' Walks the points `x_k + alpha*p_k` and returns the contrast estimate at each, `NA`
#' wherever the ESS requirement rules the point out. A point already outside `ess_bounds`
#' is skipped before `trust_function` is called at all; a point that is evaluated and
#' then fails is added to the bounds.
#'
#' The points are visited outwards from `x_k` rather than in grid order, i.e. by
#' increasing `abs(alpha_k)`. `x_k` is known to satisfy the requirement and the effective
#' sample size only falls off as the line runs away from it, so walking outwards meets
#' the breakdown on each side once and prunes the rest of that side, instead of starting
#' at one end of the line and paying for every ruled-out point before reaching the part
#' that is still usable.
#'
#' @param alpha_k Step lengths to evaluate along `p_k`.
#' @param p_k Search direction.
#' @param x_k Search starting point.
#' @param trust_function Trust function, returning a list with `f_est` and `ess`.
#' @param ess_bounds Per-parameter breakdown bounds, from \code{\link{ess_bounds_init}}.
#' @param eta_trust Smallest `ess/nSims` an iterate is allowed to have.
#' @param nSims Number of simulations the ensemble behind `trust_function` holds.
#'
#' @return A list with `z_k` (the evaluated points, one per column), `f_vals` (the
#' contrast estimate at each, `NA` where the point was ruled out) and `ess_bounds`.
trust_line_evaluate <- function(alpha_k, p_k, x_k, trust_function, ess_bounds, eta_trust, nSims){
  z_k <- matrix(data = NA, nrow = length(x_k), ncol = length(alpha_k))
  for(i in 1:length(x_k)){
    z_k[i, ] <- x_k[i] + p_k[i] * alpha_k
  }

  f_vals <- rep(NA_real_, length(alpha_k))
  for(i in order(abs(alpha_k))){
    if(!ess_bounds_feasible(ess_bounds, z_k[, i])){
      next
    }
    res <- trust_function(z_k[, i])
    passed <- (res$ess/nSims > eta_trust)
    if(passed){
      f_vals[i] <- res$f_est
    }
    ess_bounds <- ess_bounds_record(ess_bounds, z_k[, i], passed = passed)
  }

  return(list(z_k = z_k, f_vals = f_vals, ess_bounds = ess_bounds))
}

#' Per-parameter bounds on where the importance sampling breaks down
#' @description
#' A trust region step re-weights one ensemble simulated at `x_0`, so the effective
#' sample size of its estimates falls off as the evaluation point moves away from `x_0`.
#' This holds, for each parameter separately, how far the line searches have got before
#' the `ess/nSims > eta_trust` requirement failed: the smallest failing value above
#' `x_0` and the largest failing value below it. Points outside the resulting box are
#' skipped rather than evaluated, which matters because each evaluation costs a full
#' importance sampling pass over the ensemble.
#'
#' The box is a screening device, not a description of the region where the requirement
#' holds. A breakdown is a property of the whole parameter vector, and there is no way
#' to tell from one failing evaluation which parameters caused it, so the failure is
#' recorded against each of them and points that would in fact have passed can be
#' excluded. What is ruled out is therefore deliberately conservative: see
#' \code{\link{ess_bounds_record}} for the one case where it is held back.
#'
#' @param x_0 Centre of the trust region, i.e. the parameters the ensemble was simulated at.
#'
#' @return A list with `x_0`, the per-parameter `lower` and `upper` bounds (infinite until
#' a breakdown is recorded), and the per-parameter `pass_lower` and `pass_upper` extremes
#' of the values seen to satisfy the requirement.
ess_bounds_init <- function(x_0){
  # x_0 itself always passes: every weight against its own ensemble is 1, so the
  # effective sample size there is nSims exactly. It is also kept in the returned object,
  # which nothing here reads but which makes the box readable on its own -- the bounds
  # mean little without the centre they are measured from.
  return(list(x_0 = x_0,
              lower = rep(-Inf, length(x_0)),
              upper = rep(Inf, length(x_0)),
              pass_lower = x_0,
              pass_upper = x_0))
}

#' Breakdown bounds to open a step with
#' @description
#' Returns `ess_bounds` when it belongs to this `x_0`, and a fresh set otherwise. The
#' bounds are offsets measured from the centre of the trust region, so a set carried over
#' from a different centre describes a different region and is discarded rather than
#' applied to this one. That is what makes the set safe to carry across a *rejected* step,
#' where the iterate and the ensemble both stay put, but not across an accepted one.
#'
#' @param ess_bounds Bounds from a previous step, or `NULL`.
#' @param x_0 Centre of the trust region for the step about to run.
#'
#' @return An `ess_bounds` list, as \code{\link{ess_bounds_init}} returns.
ess_bounds_carried <- function(ess_bounds, x_0){
  usable <- !is.null(ess_bounds) && is.list(ess_bounds) &&
    !is.null(ess_bounds$x_0) && (length(ess_bounds$x_0) == length(x_0)) &&
    all(is.finite(ess_bounds$lower) | is.infinite(ess_bounds$lower)) &&
    isTRUE(all.equal(ess_bounds$x_0, x_0))
  if(usable){
    return(ess_bounds)
  }
  return(ess_bounds_init(x_0))
}

#' Record the outcome of one ESS requirement check
#' @description
#' A failing point tightens the bound of each parameter on the side it moved to, but
#' only past values of that parameter no evaluation has yet passed at. Without that
#' restriction a failure driven by one parameter drags the bounds of the others in to
#' wherever they happened to sit, which ratchets the box shut around the search: the
#' line searches move a parameter out to near its own limit, a failure elsewhere then
#' records that near-limit value as the limit, and the parameter can never return to it.
#'
#' Holding the bounds outside every value that has passed also keeps every point the
#' search has accepted inside the box, which is what lets a line search fall back on its
#' starting point when the rest of the line is ruled out.
#'
#' @param ess_bounds Bounds to record into, from \code{\link{ess_bounds_init}}.
#' @param x Parameter vector that was evaluated.
#' @param passed TRUE when `x` satisfied `ess/nSims > eta_trust`.
#'
#' @return `ess_bounds` updated with the outcome.
ess_bounds_record <- function(ess_bounds, x, passed){
  if(passed){
    ess_bounds$pass_lower <- pmin(ess_bounds$pass_lower, x)
    ess_bounds$pass_upper <- pmax(ess_bounds$pass_upper, x)
    return(ess_bounds)
  }

  # Parameters the failing point left at, or inside, their passing range are not
  # evidence of anything and keep the bounds they have. This also covers the parameters
  # the point never moved off x_0.
  above <- (x > ess_bounds$pass_upper)
  below <- (x < ess_bounds$pass_lower)
  ess_bounds$upper[above] <- pmin(ess_bounds$upper[above], x[above])
  ess_bounds$lower[below] <- pmax(ess_bounds$lower[below], x[below])
  return(ess_bounds)
}

#' Check a point against the recorded ESS breakdown bounds
#' @description
#' The bounds are values at which the requirement was seen to fail, so a point sitting
#' exactly on one is itself ruled out.
#'
#' @param ess_bounds Bounds to check against, from \code{\link{ess_bounds_init}}.
#' @param x Parameter vector to check.
#'
#' @return TRUE when `x` is not already known to break the ESS requirement.
ess_bounds_feasible <- function(ess_bounds, x){
  return(all((x > ess_bounds$lower) & (x < ess_bounds$upper)))
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
  const <- t(x_k-x_0)%*%(x_k-x_0) - delta^2
  # The discriminant is 4*a*(delta^2 - d^2), where d is the perpendicular distance from
  # x_0 to the line, so it cannot really be negative while x_k is inside the region. It
  # can still come out a hair below zero by rounding when x_k sits essentially exactly on
  # the boundary and p_k is near-tangent, and the NaN then propagates into seq() as an
  # error that aborts the whole fit, so it is clamped rather than left to escape.
  determ <- sqrt(max(b^2-4*a*const, 0))
  up_lim <- (determ - b)/(2*a)
  low_lim <- (-1)*(determ + b)/(2*a)
  return(c(min(c(low_lim, up_lim)), max(c(low_lim, up_lim))))
}
