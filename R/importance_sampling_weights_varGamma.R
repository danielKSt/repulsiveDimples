#' Calculate importance sampling weights for a Variance gamma process
#' @description
#' Calculate importance sampling weights for a Variance gamma process, for moving
#' `(kappa, scale, mu, nu)` away from the `(kappa_0, scale_0, mu_0, nu_0)` the ensemble
#' was simulated with. `nu_0` defaults to `nu`, in which case the shape is not reweighted
#' and `nu` has to be the shape the ensemble was simulated with.
#'
#' @param kappa Kappa parameter
#' @param mu mu parameter
#' @param scale scale parameter, \eqn{\sigma} in \code{\link{varGamma_daughter_log_density}}
#' @param nu Shape parameter of the kernel, greater than -1
#' @param kappa_0 Kappa parameter the ensemble was simulated with
#' @param mu_0 mu parameter the ensemble was simulated with
#' @param scale_0 scale parameter the ensemble was simulated with
#' @param nu_0 Shape parameter the ensemble was simulated with, greater than -1. Defaults
#' to `nu`.
#' @param log_f_kappa_0 Parent null density
#' @param log_fCond_theta_0 Daughter null density
#' @param patternSim List of point patterns to use for estimation
#' @param parallel Set to TRUE to use parallel computing
#' @param daughter_kernel_cache Optional cache from a previous call on the same ensemble,
#' as returned in this call's output: a list with `scale`, `nu`, per-pattern
#' `kernel_sums`, and `window_grid` (see \code{\link{varGamma_window_grid}}). The kernel
#' sums, which depend on `scale` and `nu` but not on `mu`, are only recomputed when `scale`
#' or `nu` differs from the cached value. The window grid does not depend on the
#' parameters at all, so it is kept across changes of `scale` and `nu`, and only extended
#' when they reach outside what it covers. A grid built for a different ensemble is
#' detected and rebuilt. Pass `NULL` (the default) to start from scratch.
#'
#' @return A list with `w_is` (the importance sampling weights) and `daughter_kernel_cache`
#' (to be passed back into the next call). When the target and simulation parameters
#' coincide the weights are all exactly 1 and `daughter_kernel_cache` is returned unchanged.
#'
#' @export
importance_sampling_weights_varGamma <- function(kappa, mu, scale, nu,
                                                 kappa_0 = NULL, mu_0 = NULL, scale_0 = NULL, nu_0 = nu,
                                                 patternSim, log_f_kappa_0 = NULL, log_fCond_theta_0 = NULL,
                                                 parallel = FALSE, daughter_kernel_cache = NULL){

  if(isTRUE(kappa == kappa_0) && isTRUE(scale == scale_0) && isTRUE(mu == mu_0) && isTRUE(nu == nu_0)){
    return(list(w_is = rep(1, length(patternSim)),
                daughter_kernel_cache = daughter_kernel_cache)) # All parameters being equal gives w_is = 1
  }
  lapply_fun <- if(parallel) parallel::mclapply else lapply

  # Parent density
  B_area   <- vapply(patternSim, function(p) p$B_area, numeric(1))
  n_parent <- vapply(patternSim, function(p) nrow(p$parent), numeric(1))

  if(is.null(log_f_kappa_0)){
    log_f_kappa_0 <- B_area*(1 - kappa_0) + n_parent*log(kappa_0)
  }
  log_f_kappa <- B_area*(1 - kappa) + n_parent*log(kappa)

  # Daughters density
  # The kernel sums depend on scale and nu only, so they are computed once per (scale, nu)
  # and reused. The window grid behind their window probabilities is kept even when
  # (scale, nu) moves, since it depends on the ensemble alone.
  window_grid <- daughter_kernel_cache$window_grid
  if(is.null(daughter_kernel_cache) || !isTRUE(all.equal(daughter_kernel_cache$scale, scale)) ||
     !isTRUE(all.equal(daughter_kernel_cache$nu, nu))){
    sums <- varGamma_ensemble_kernel_sums(patternSim = patternSim, scale = scale, nu = nu,
                                          window_grid = window_grid, lapply_fun = lapply_fun)
    window_grid <- sums$window_grid
    daughter_kernel_cache <- list(scale = scale, nu = nu, kernel_sums = sums$kernel_sums,
                                  window_grid = window_grid)
  }

  log_fCond_theta <- vapply(daughter_kernel_cache$kernel_sums,
                            daughter_log_density_from_sums, numeric(1), mu = mu)

  if(is.null(log_fCond_theta_0)){
    if(isTRUE(all.equal(scale_0, scale)) && isTRUE(all.equal(nu_0, nu))){
      # Baseline and target share the kernel, so the kernel sums just computed are the
      # ones the baseline density needs as well, and only the mu-dependent combine differs.
      log_fCond_theta_0 <- vapply(daughter_kernel_cache$kernel_sums,
                                  daughter_log_density_from_sums, numeric(1), mu = mu_0)
    } else {
      # The same window grid as the target, so that both sides of the weight carry the
      # same approximation of the window probabilities.
      sums_0 <- varGamma_ensemble_kernel_sums(patternSim = patternSim, scale = scale_0, nu = nu_0,
                                              window_grid = daughter_kernel_cache$window_grid,
                                              lapply_fun = lapply_fun)
      daughter_kernel_cache$window_grid <- sums_0$window_grid
      log_fCond_theta_0 <- vapply(sums_0$kernel_sums, daughter_log_density_from_sums,
                                  numeric(1), mu = mu_0)
    }
  }

  # Return final importance sampling weights
  w_is <- exp(log_f_kappa + log_fCond_theta - log_f_kappa_0 - log_fCond_theta_0)
  return(list(w_is = w_is, daughter_kernel_cache = daughter_kernel_cache))
}

#' Log-density of the daughter pattern given parents for a Variance gamma process
#'
#' @description
#' Computes the log-density of a daughter point pattern \eqn{X} given a parent point
#' pattern \eqn{C}, under a Variance gamma cluster process. Conditional on \eqn{C}, the
#' daughters form a Poisson point process on \eqn{W} with intensity
#' \eqn{\mu \sum_{c_j \in C} k(x - c_j)}, where
#' \deqn{k(u) = \frac{(\|u\|/\sigma)^{\nu} K_{\nu}(\|u\|/\sigma)}{\pi 2^{\nu+1} \sigma^2 \Gamma(\nu+1)}}
#' is the Variance gamma kernel of Jalilian, Guan and Waagepetersen (2013), with
#' \eqn{K_\nu} the modified Bessel function of the second kind. At \eqn{\nu = 1/2} this is
#' the exponential kernel \eqn{k(u) = \exp(-\|u\|/\sigma)/(2\pi\sigma^2)}. The log-density
#' with respect to the unit rate Poisson process is
#' \deqn{
#'   \log f(X \mid C) = |W| - \mu \sum_{c_j \in C} p_W(c_j)
#'   + \sum_{x_i \in X} \log\left(\mu \sum_{c_j \in C} k(x_i - c_j)\right)
#' }
#' where \eqn{p_W(c_j) = \int_W k(x - c_j)\,\mathrm{d}x} is the probability that a
#' daughter of \eqn{c_j} falls inside \eqn{W}, see \code{\link{varGamma_daughter_kernel_sums}}
#' for how it is computed.
#'
#' @param mu       Positive scalar. Mean number of daughters per parent.
#' @param scale    Positive scalar. The scale \eqn{\sigma} of the kernel, `scale` in
#' \code{\link{rVarGamma_matern_thinned}}.
#' @param nu       Scalar greater than -1. The shape \eqn{\nu} of the kernel.
#' @param parent   A \code{data.frame} with columns \code{x} and \code{y} giving parent locations,
#' typically restricted to an enlarged window \eqn{B \supset W}.
#' @param daughter A \code{data.frame} with columns \code{x} and \code{y} giving daughter locations,
#' restricted to the observation window \eqn{W}.
#' @param xlim     Numeric vector of length 2. The x-range of the observation window \eqn{W}.
#' @param ylim     Numeric vector of length 2. The y-range of the observation window \eqn{W}.
#' @param pattern  A list containing elements "parent", "daughter", "xlim_unthinned", and "ylim_unthinned".
#'
#' @return A scalar giving \eqn{\log f(X \mid C)}.
#'
#' @export
varGamma_daughter_log_density <- function(pattern = NULL, mu, scale, nu,
                                          parent = NULL, daughter = NULL,
                                          xlim = NULL, ylim = NULL){
  kernel_sums <- varGamma_daughter_kernel_sums(pattern = pattern, scale = scale, nu = nu,
                                               parent = parent, daughter = daughter,
                                               xlim = xlim, ylim = ylim)
  return(daughter_log_density_from_sums(mu = mu, kernel_sums = kernel_sums))
}

#' Kernel sums for the Variance gamma daughter log-density (the scale and nu dependent part)
#'
#' @description
#' Computes the parts of \code{\link{varGamma_daughter_log_density}} that depend on the
#' kernel but not on \eqn{\mu}: \eqn{\sum_j p_W(c_j)} and
#' \eqn{\sum_i \log \sum_j k(x_i-c_j)}. They are combined with \eqn{\mu} by
#' \code{\link{daughter_log_density_from_sums}}.
#'
#' The 2D Variance gamma kernel does not factorise over the two coordinates, so unlike
#' for the Thomas process \eqn{p_W(c)} has no closed form on a rectangle. It is computed
#' from the normal variance mixture the kernel is simulated from: a daughter is displaced
#' from its parent by \eqn{\sqrt{S} Z}, with \eqn{Z} standard normal in 2D and
#' \eqn{S = 2\sigma^2 G}, \eqn{G \sim \mathrm{Gamma}(\nu + 1, 1)}. Given \eqn{S = s} the
#' probability \eqn{I(c; s)} that the daughter lands in \eqn{W} factorises into normal CDF
#' differences, as for the Thomas process. Only the sum over the parents is needed, and
#' \eqn{J(s) = \sum_j I(c_j; s)} does not depend on \eqn{(\sigma, \nu)}, so the sum is
#' computed from a precomputed table of \eqn{J}, see \code{\link{varGamma_window_grid}}
#' for how and how accurately. This function builds that table for its one pattern on
#' the fly; \code{\link{importance_sampling_weights_varGamma}} builds it once per
#' ensemble instead and reuses it for every \eqn{(\sigma, \nu)}.
#'
#' @param scale    Positive scalar. The scale \eqn{\sigma} of the kernel.
#' @param nu       Scalar greater than -1. The shape \eqn{\nu} of the kernel.
#' @param parent   A \code{data.frame} with columns \code{x} and \code{y} giving parent locations.
#' @param daughter A \code{data.frame} with columns \code{x} and \code{y} giving daughter locations.
#' @param xlim     Numeric vector of length 2. The x-range of the observation window \eqn{W}.
#' @param ylim     Numeric vector of length 2. The y-range of the observation window \eqn{W}.
#' @param pattern  A list containing elements "parent", "daughter", "xlim_unthinned", and "ylim_unthinned".
#'
#' @return A list with \code{area_W}, \code{S_pW} (\eqn{\sum_j p_W(c_j)}),
#' \code{sumLogK} (\eqn{\sum_i \log \sum_j k(x_i-c_j)}), and \code{nDaughter}.
#'
#' @export
varGamma_daughter_kernel_sums <- function(pattern = NULL, scale, nu,
                                          parent = NULL, daughter = NULL,
                                          xlim = NULL, ylim = NULL){
  if(is.null(pattern)){
    pattern <- list(parent = parent, daughter = daughter,
                    xlim_unthinned = xlim, ylim_unthinned = ylim)
  }
  sums <- varGamma_ensemble_kernel_sums(patternSim = list(pattern), scale = scale, nu = nu)
  return(sums$kernel_sums[[1]])
}

#' Window grid for the Variance gamma window probabilities
#'
#' @description
#' Tabulates, for every pattern of an ensemble, the part of the window probabilities
#' \eqn{\sum_{c \in C} p_W(c)} that does not depend on the kernel parameters, so that
#' \code{\link{varGamma_window_sums}} can then give the sums at any \eqn{(\sigma, \nu)}
#' with one matrix-vector product.
#'
#' With \eqn{S = 2\sigma^2 G}, \eqn{G \sim \mathrm{Gamma}(\nu + 1, 1)}, the mixing variance of
#' the kernel (see \code{\link{varGamma_daughter_kernel_sums}}), and \eqn{I(c; s)} the
#' probability that a normal displacement with covariance \eqn{sI} from \eqn{c} lands in
#' \eqn{W},
#' \deqn{\sum_{c} p_W(c) = \int_0^\infty \gamma(s)\, J(s)\,\mathrm{d}s, \qquad J(s) = \sum_{c} I(c; s),}
#' with \eqn{\gamma} the Gamma\eqn{(\nu+1, 1/(2\sigma^2))} density. \eqn{J} depends on the
#' parents and the window only. As \eqn{s \to 0}, \eqn{J(s)} tends to \eqn{J_0}, the number
#' of parents inside \eqn{W}, so after substituting \eqn{t = \log s},
#' \deqn{\sum_{c} p_W(c) = J_0 + \int_{-\infty}^{\infty} g(t)\,\bigl(J(e^t) - J_0\bigr)\,\mathrm{d}t,
#'   \qquad g(t) = \frac{b^a}{\Gamma(a)} e^{a t - b e^t},}
#' with \eqn{a = \nu + 1} and \eqn{b = 1/(2\sigma^2)}. Subtracting \eqn{J_0} makes the integrand
#' vanish quickly at both ends whatever \eqn{\nu}: on the left because \eqn{J(e^t) - J_0}
#' does, like \eqn{\exp(-d^2 e^{-t}/2)} for a parent at distance \eqn{d} from an edge, even
#' where \eqn{g} itself, for \eqn{\nu} near -1, decays only slowly; on the right because
#' \eqn{g} does. The integrand is also smooth in \eqn{t}, so the trapezoidal rule on a
#' uniform grid in \eqn{t} converges geometrically as the spacing shrinks. The grid stores
#' \eqn{D_k = J(e^{t_k}) - J_0} per pattern, and the sums are
#' \eqn{J_0 + h \sum_k g(t_k) D_k}.
#'
#' The grid starts at \eqn{t = \log(d^2/100)}, with \eqn{d} the smallest distance from any
#' parent of the ensemble to an edge line of its window (bounded below by \eqn{10^{-12}}
#' times the window side), below which \eqn{D_k} is negligible. It ends where the upper
#' \eqn{10^{-20}} tail of \eqn{S} begins, plus a margin of 2. The spacing starts at `h` and is
#' halved until \eqn{(\nu + 1) h^2 \le 3/8}, since \eqn{g} narrows like \eqn{1/\sqrt{\nu+1}}.
#' The resulting relative error of the sums is about \eqn{10^{-11}} for
#' \eqn{-0.95 \le \nu \le 20}. Since the left end and the table do not depend on
#' \eqn{(\sigma, \nu)}, a grid built at one \eqn{(\sigma, \nu)} is extended to the right, or
#' refined, when later used at another that needs it.
#'
#' @param patternSim List of point patterns, each with `parent`, `xlim_unthinned` and
#' `ylim_unthinned`.
#' @param scale Positive scalar. The scale \eqn{\sigma} the grid has to cover at first.
#' @param nu Scalar greater than -1. The shape \eqn{\nu} the grid has to cover at first.
#' @param h Initial spacing of the grid in \eqn{t = \log s}.
#' @param lapply_fun The `lapply` to tabulate the patterns with, for instance
#' \code{\link[parallel]{mclapply}}.
#'
#' @return A list with `t` (the grid), `h` (its spacing), `J0` (\eqn{J_0} per pattern), `D`
#' (a matrix with one row of \eqn{D_k} per pattern and one column per grid point), and
#' `fingerprint` (identifying the ensemble the grid was built for).
#'
varGamma_window_grid <- function(patternSim, scale, nu, h = 0.25, lapply_fun = lapply){
  if(!(nu > -1)){
    stop("nu must be greater than -1")
  }
  side <- max(vapply(patternSim, function(p) max(diff(p$xlim_unthinned), diff(p$ylim_unthinned)),
                     numeric(1)))
  d_min <- min(vapply(patternSim, function(p){
    if(nrow(p$parent) == 0) return(Inf)
    xl <- p$xlim_unthinned; yl <- p$ylim_unthinned
    min(abs(p$parent$x - xl[1]), abs(p$parent$x - xl[2]),
        abs(p$parent$y - yl[1]), abs(p$parent$y - yl[2]))
  }, numeric(1)))
  # A parent closer to an edge than this is a parent on the edge to within the precision of
  # its coordinates, and an ensemble without parents has nothing to tabulate.
  d_min <- if(is.finite(d_min)) max(d_min, 1e-12*side) else side

  h <- varGamma_window_grid_spacing(h, nu)
  t_min <- log(d_min^2/100)
  t_max <- max(varGamma_window_grid_t_max(scale, nu) + 2, t_min + h)
  t <- t_min + h*(0:ceiling((t_max - t_min)/h))

  rows <- varGamma_window_grid_rows(patternSim, t, lapply_fun)
  return(list(t = t, h = h, J0 = rows$J0, D = rows$D,
              fingerprint = varGamma_window_grid_fingerprint(patternSim)))
}

#' Window probability sums from a Variance gamma window grid
#'
#' @description
#' Evaluates \eqn{\sum_{c \in C} p_W(c)} for every pattern a window grid was built for, at
#' the kernel parameters \eqn{(\sigma, \nu)}, by the trapezoidal rule described in
#' \code{\link{varGamma_window_grid}}. The grid has to cover \eqn{(\sigma, \nu)}; one that
#' was built at, or extended to, these parameters does.
#'
#' @param window_grid A grid from \code{\link{varGamma_window_grid}}.
#' @param scale Positive scalar. The scale \eqn{\sigma} of the kernel.
#' @param nu Scalar greater than -1. The shape \eqn{\nu} of the kernel.
#'
#' @return A vector with \eqn{\sum_{c \in C} p_W(c)} for each pattern.
#'
varGamma_window_sums <- function(window_grid, scale, nu){
  a <- nu + 1
  b <- 1/(2*scale^2)
  t <- window_grid$t
  w <- window_grid$h * exp(a*log(b) - lgamma(a) + a*t - b*exp(t))
  return(window_grid$J0 + drop(window_grid$D %*% w))
}

#' Make a Variance gamma window grid cover the kernel parameters
#'
#' @description
#' Returns a window grid for `patternSim` that covers `(scale, nu)`: `window_grid` itself
#' when it already does, `window_grid` refined or extended to the right when it does not,
#' and a fresh grid from \code{\link{varGamma_window_grid}} when there is none or it was
#' built for another ensemble.
#'
#' @inheritParams varGamma_window_grid
#' @param window_grid A grid from \code{\link{varGamma_window_grid}}, or `NULL`.
#'
#' @return A window grid, see \code{\link{varGamma_window_grid}}.
varGamma_window_grid_cover <- function(window_grid, patternSim, scale, nu, lapply_fun = lapply){
  if(!(nu > -1)){
    stop("nu must be greater than -1")
  }
  if(is.null(window_grid) ||
     !identical(window_grid$fingerprint, varGamma_window_grid_fingerprint(patternSim))){
    return(varGamma_window_grid(patternSim = patternSim, scale = scale, nu = nu,
                                lapply_fun = lapply_fun))
  }
  # Refine, by tabulating the midpoints, until the spacing resolves g at this nu.
  while(varGamma_window_grid_spacing(window_grid$h, nu) < window_grid$h){
    t_mid <- window_grid$t[-length(window_grid$t)] + window_grid$h/2
    rows <- varGamma_window_grid_rows(patternSim, t_mid, lapply_fun)
    t_all <- c(window_grid$t, t_mid)
    ord <- order(t_all)
    window_grid$t <- t_all[ord]
    window_grid$D <- cbind(window_grid$D, rows$D)[, ord, drop = FALSE]
    window_grid$h <- window_grid$h/2
  }
  # Extend to the right until the grid reaches past the bulk of the mixing density.
  t_last <- window_grid$t[length(window_grid$t)]
  t_need <- varGamma_window_grid_t_max(scale, nu)
  if(t_last < t_need){
    t_new <- t_last + window_grid$h*seq_len(ceiling((t_need + 2 - t_last)/window_grid$h))
    rows <- varGamma_window_grid_rows(patternSim, t_new, lapply_fun)
    window_grid$t <- c(window_grid$t, t_new)
    window_grid$D <- cbind(window_grid$D, rows$D)
  }
  return(window_grid)
}

#' Kernel sums for every pattern of a Variance gamma ensemble
#'
#' @description
#' \code{\link{varGamma_daughter_kernel_sums}} for a whole ensemble, with the window
#' probabilities taken from a window grid shared by the ensemble (see
#' \code{\link{varGamma_window_grid}}), which is made to cover `(scale, nu)` first and
#' returned for reuse.
#'
#' @inheritParams varGamma_window_grid_cover
#'
#' @return A list with `kernel_sums` (one list per pattern, as returned by
#' \code{\link{varGamma_daughter_kernel_sums}}) and `window_grid`.
varGamma_ensemble_kernel_sums <- function(patternSim, scale, nu, window_grid = NULL, lapply_fun = lapply){
  window_grid <- varGamma_window_grid_cover(window_grid = window_grid, patternSim = patternSim,
                                            scale = scale, nu = nu, lapply_fun = lapply_fun)
  S_pW <- varGamma_window_sums(window_grid, scale = scale, nu = nu)
  sumLogK <- lapply_fun(patternSim, function(p){
    varGamma_log_kernel_sum(parent = p$parent, daughter = p$daughter, scale = scale, nu = nu)
  })
  kernel_sums <- lapply(seq_along(patternSim), function(i){
    p <- patternSim[[i]]
    list(area_W = diff(p$xlim_unthinned)*diff(p$ylim_unthinned), S_pW = S_pW[i],
         sumLogK = sumLogK[[i]], nDaughter = nrow(p$daughter))
  })
  return(list(kernel_sums = kernel_sums, window_grid = window_grid))
}

# \sum_i \log \sum_j k(x_i - c_j) for the Variance gamma kernel.
varGamma_log_kernel_sum <- function(parent, daughter, scale, nu){
  if(nrow(daughter) == 0){
    return(0)
  }
  dx <- outer(daughter$x, parent$x, "-")
  dy <- outer(daughter$y, parent$y, "-")
  z <- sqrt(dx*dx + dy*dy)/scale
  # The kernel is z^nu K_nu(z) up to a constant factor, which comes out of the row sums
  # as a single term per daughter.
  E <- z^nu * besselK(z, nu)
  sum(log(rowSums(E))) - nrow(daughter)*log(pi * 2^(nu + 1) * scale^2 * gamma(nu + 1))
}

# J_0 and the rows D_k = J(exp(t_k)) - J_0 of a window grid, one row per pattern.
#' @importFrom stats pnorm
varGamma_window_grid_rows <- function(patternSim, t, lapply_fun = lapply){
  sd <- exp(t/2)
  rows <- lapply_fun(patternSim, function(p){
    if(nrow(p$parent) == 0){
      return(list(J0 = 0, D = numeric(length(t))))
    }
    xl <- p$xlim_unthinned; yl <- p$ylim_unthinned
    cx <- p$parent$x; cy <- p$parent$y
    # one row per parent, one column per grid point
    Px <- stats::pnorm(outer(xl[2] - cx, sd, "/")) - stats::pnorm(outer(xl[1] - cx, sd, "/"))
    Py <- stats::pnorm(outer(yl[2] - cy, sd, "/")) - stats::pnorm(outer(yl[1] - cy, sd, "/"))
    J0 <- sum(cx > xl[1] & cx < xl[2] & cy > yl[1] & cy < yl[2])
    list(J0 = J0, D = colSums(Px*Py) - J0)
  })
  return(list(J0 = vapply(rows, function(r) r$J0, numeric(1)),
              D = matrix(unlist(lapply(rows, function(r) r$D)), nrow = length(rows),
                         byrow = TRUE)))
}

# Largest grid spacing, no larger than h, that resolves the mixing density at nu.
varGamma_window_grid_spacing <- function(h, nu){
  while((nu + 1)*h^2 > 3/8){
    h <- h/2
  }
  h
}

# log s beyond which the mixing density of s has only its upper 1e-20 tail left.
#' @importFrom stats qgamma
varGamma_window_grid_t_max <- function(scale, nu){
  log(2*scale^2*stats::qgamma(1e-20, shape = nu + 1, lower.tail = FALSE))
}

# Identifies the ensemble a window grid belongs to: the parent count and coordinate sums
# of every pattern, and its window.
varGamma_window_grid_fingerprint <- function(patternSim){
  vapply(patternSim, function(p){
    c(nrow(p$parent), if(nrow(p$parent) > 0) c(sum(p$parent$x), sum(p$parent$y)) else c(0, 0),
      p$xlim_unthinned, p$ylim_unthinned)
  }, numeric(7))
}
