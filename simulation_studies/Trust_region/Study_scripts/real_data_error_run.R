# How much of the error in a fitted parameter is the optimizer, and how much is estimating
# K and rho in the first place?
#
# Every trust region study so far has fitted against study1's K_hat and rho_hat, which are
# averaged over nTrue = 1e6 realisations and so are effectively the exact population
# targets. Those studies therefore measure optimizer error with *no* estimation error in
# the estimand at all. On real data the targets are estimated from a handful of snapshots,
# and the claim that the optimizer's remaining error is immaterial in practice rests on
# that estimation error being the larger of the two. This measures it.
#
# Design: for each replicate, simulate `m` fresh snapshots at the study's own window size,
# estimate rho and K from them exactly as an application would (intensity_est and
# K_est.unions, the union-of-windows estimators), derive the starting parameters from that
# same estimated K with thomas.estK, and fit. The whole chain is driven by the m snapshots,
# so the error that comes back is what an analyst with m snapshots would actually get.
# Scoring is against par_thomas, so it is directly comparable with every other study.
#
# The comparison number is the study1_quad error at the same nSims, which is the m = Inf
# case -- the same optimizer fitting the exact targets. If total error at realistic m is
# far above it, the optimizer is not the binding constraint on real data and the
# recommendation to prefer the cheaper method follows. If it is close, it does not.
#
# m = 1 is the case that matters most: a single snapshot is the commonest real-data
# situation, and the one where estimating the targets costs the most. It relies on
# K_est.unions handling one snapshot, which it did not before 38c3e82 -- two of its loops
# counted backwards, duplicating the pattern and giving K(1.0) = 9.148 where the same data
# at m = 5 gives 3.824. Checked after the fix: at m = 1 it now agrees exactly with a
# direct Kest on that snapshot.

library(spatstat)
library(parallel)
library(mcprogress)
library(repulsiveDimples)

resFolder <- "simulation_studies/Trust_region/Results/"
dir.create(resFolder, showWarnings = FALSE, recursive = TRUE)
load(file = "simulation_studies/Trust_region/Data/study1.RDa")

r_vec   <- seq(from = 0, to = 3, by = 0.05)   # the grid study1's K_hat is on
spacing <- 5                                  # as in Data/study1_data.R
nSims   <- 5000        # optimizer ensemble; study1_quad median error at this level is 0.056
mValues <- c(1, 2, 5, 10, 25, 50)
nReps   <- 25
nCores  <- 32

# study1_quad medians, for reading the results against: the same optimizer on exact targets.
optimizerOnly <- c("500" = 0.097, "1000" = 0.076, "5000" = 0.046, "10000" = 0.044)

truth <- log(c(par_thomas$kappa, par_thomas$omega, par_thomas$mu))
dr    <- r_vec[2] - r_vec[1]

# One replicate: m snapshots in, one fitted parameter vector out.
one_rep <- function(i, m, seedBase){
  options(mc.cores = 1)
  set.seed(seedBase + i)

  sims <- lapply(seq_len(m), function(j)
    rThomas_matern_thinned(kappa = par_thomas$kappa, scale = par_thomas$omega,
                           mu = par_thomas$mu, repulsionRange = par_thomas$rRange,
                           xlims = c(0, sidelength), ylims = c(0, sidelength),
                           saveparents = TRUE))
  points_input <- lapply(sims, function(a) if(nrow(a$thinned) == 0) 0 else a$thinned)

  rho_m <- intensity_est(points_input = points_input, l = sidelength, timescale = 1)
  K_m   <- K_est.unions(points_input = points_input, l = sidelength, spacing = spacing,
                        r_vec = r_vec, timescale = 1)

  # How far the estimated targets themselves are from the population ones. The K distance
  # is the same functional the contrast minimises, so it is on the scale the fit cares
  # about rather than an arbitrary norm.
  rho_err <- abs(rho_m - rho_hat)/rho_hat
  K_err   <- sqrt(sum((K_m$border^(1/4) - K_hat$border^(1/4))^2, na.rm = TRUE)*dr)

  # The starting point is derived from this replicate's own K, as it would be in practice,
  # so its error propagates into the fit rather than being handed in for free.
  st <- try(thomas.estK(X = K_m, rmin = 2*par_thomas$rRange, rmax = 3)$par, silent = TRUE)
  if(inherits(st, "try-error") || any(!is.finite(st)) || any(st <= 0)){
    return(list(m = m, rep = i, par = rep(NA_real_, 3), start = rep(NA_real_, 3),
                rho_err = rho_err, K_err = K_err, secs = NA_real_, note = "thomas.estK failed"))
  }
  lambda <- -log(1 - pi*rho_m*par_thomas$rRange^2)/(pi*par_thomas$rRange^2)
  start_params <- log(c(st[1], st[2], lambda/st[1]))
  names(start_params) <- c("kappa", "omega", "mu")
  if(any(!is.finite(start_params))){
    return(list(m = m, rep = i, par = rep(NA_real_, 3), start = rep(NA_real_, 3),
                rho_err = rho_err, K_err = K_err, secs = NA_real_, note = "bad start"))
  }

  t0 <- proc.time()
  fit <- try(min_contrast_trust_region(
    params = start_params, parFreeIndex = c(1, 2, 3), repRange = par_thomas$rRange,
    rho_hat = rho_m, K_hat = K_m, xlims = c(0, 3), ylims = c(0, 3), nSims = nSims,
    deltaInit = 0.2, eta = 0.05, deltaMax = 1.0, deltaMin = 0.0001, wq = c(1000, 1/4),
    normalized = TRUE, eta_trust = 0.6, eta_converged = 0.999, validate.converged = FALSE,
    method.main = "quadratic", method.prefit = "quadratic",
    maxPrefitLoops = 8, max.iter.prefit = 10, tol = 10^-8, max.iter = 30,
    printProgress = FALSE), silent = TRUE)
  el <- proc.time() - t0
  if(inherits(fit, "try-error")){
    return(list(m = m, rep = i, par = rep(NA_real_, 3), start = as.numeric(start_params),
                rho_err = rho_err, K_err = K_err, secs = NA_real_, note = "fit failed"))
  }
  list(m = m, rep = i, par = as.numeric(fit$params[nrow(fit$params), ]),
       start = as.numeric(start_params), rho_err = rho_err, K_err = K_err,
       secs = unname(el[["elapsed"]]), note = "ok")
}

res <- list()
# Seeds are keyed on m rather than on its position in mValues, so editing the list does
# not reseed the levels that stay in it. With nReps <= 999 no two levels overlap.
for(mi in seq_along(mValues)){
  m <- mValues[mi]
  message(sprintf("\n=== m = %d snapshots, %d replicates ===", m, nReps))
  t0 <- Sys.time()
  res[[as.character(m)]] <- mcprogress::pmclapply(seq_len(nReps), one_rep, m = m,
                                                  seedBase = 80000 + 1000*m,
                                                  mc.cores = nCores, mc.preschedule = FALSE)
  message(sprintf("m = %3d done in %.1f min", m,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  save(res, mValues, nReps, nSims, r_vec, spacing, truth, par_thomas, rho_hat, K_hat,
       optimizerOnly, file = paste0(resFolder, "study1_real_data_error.RDa"))
}

# Results: ----
L2 <- function(v) sqrt(sum((v - truth)^2))
base <- optimizerOnly[[as.character(nSims)]]
cat(sprintf("\n\nOptimizer-only error at nSims = %d (study1_quad, exact targets): %.4f\n\n", nSims, base))
cat("   m   n   rho err   K err   start err   FINAL err   implied estimation part\n")
for(m in mValues){
  r  <- Filter(function(x) x$note == "ok", res[[as.character(m)]])
  if(length(r) == 0){ cat(sprintf("%4d   0   (every replicate failed)\n", m)); next }
  ef <- vapply(r, function(x) L2(x$par), numeric(1))
  es <- vapply(r, function(x) L2(x$start), numeric(1))
  # Treating the two sources as roughly independent, which is a simplification: the
  # estimated targets also move where the optimum *is*, they do not merely add noise.
  implied <- if(median(ef) > base) sqrt(median(ef)^2 - base^2) else NA_real_
  cat(sprintf("%4d  %2d    %.3f   %.4f     %.4f      %.4f      %s\n", m, length(r),
      median(vapply(r, function(x) x$rho_err, numeric(1))),
      median(vapply(r, function(x) x$K_err,  numeric(1))),
      median(es), median(ef),
      if(is.na(implied)) "below optimizer error" else sprintf("%.4f  (%.1fx optimizer)", implied, implied/base)))
}
failed <- sum(vapply(unlist(res, recursive = FALSE), function(x) x$note != "ok", logical(1)))
if(failed > 0) cat(sprintf("\n%d replicate(s) did not produce a fit; see the `note` field.\n", failed))
cat("\nRead the FINAL column against the optimizer-only number at the top. If it is far\n")
cat("larger at a realistic m, the optimizer is not what limits accuracy on real data.\n")
