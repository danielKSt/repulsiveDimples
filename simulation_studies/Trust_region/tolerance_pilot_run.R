# Paired pilot: does a smaller tolPrefitLoops give a better fit?
#
# Motivation, measured on study1_quad.RDa (2026-09-21). The main loop is close to
# useless: across both methods and all four ensemble sizes it improved the fit in only
# 40-57% of fits, i.e. a coin flip, and the final error correlates 0.53-0.93 with where
# pre-fitting left off. So the fit is decided in the pre-fitting blocks. Within that same
# run, fits that happened to use more pre-fit cycles ended better -- Spearman rho -0.23 to
# -0.32, p <= 0.001 at every level -- and at nSims = 10000 the 19 fits that ran 5 cycles
# reached a median error of 0.022 against 0.044 at 4 cycles. That is observational: a fit
# runs more cycles because its cycles kept moving. This pilot randomises nothing either,
# but it does hold the seed fixed and vary only the tolerance, which the earlier
# comparison could not.
#
# Design: same seeds in both arms, so fit i sees the same starting ensemble whichever
# tolerance it is run under and the comparison is within-fit. nSims = 10000 because that
# is where the plateau is clearest and where the conjugate main loop does still work.
#
# Both methods are included because the mechanism is not specific to either: if the fit is
# decided by pre-fitting then the tolerance moves both, and this is not a fix that closes
# the gap between the two methods. Arms run quadratic first and each arm is saved as it
# finishes, so the quadratic answer is in hand before the conjugate arms are done and the
# run can be stopped early without losing anything.
#
# Cost: at the observed 1854 s per quadratic fit at nSims = 10000, and more for the
# smaller tolerance since it buys more cycles, expect roughly 2-2.5 h for the two
# quadratic arms on 32 cores and a similar amount again for the conjugate pair.

library(spatstat)
library(parallel)
library(mcprogress)
library(repulsiveDimples)

source("simulation_studies/Trust_region/study_run_helpers.R")
dir.create(resFolder, showWarnings = FALSE, recursive = TRUE)

load(file = "simulation_studies/Trust_region/Data/study1.RDa")
start_params <- study_start_params(K_hat_unions = K_hat_unions, par_thomas = par_thomas,
                                   rho_hat = rho_hat)
truth <- log(c(par_thomas$kappa, par_thomas$omega, par_thomas$mu))
rm(K_hat_unions, nTrue, sidelength, nUnions)

nRuns  <- 40
nSims  <- 10000
seed   <- 70000          # shared by every arm, so fit i is paired across all four
tols   <- c(0.05, 0.02)  # current default against the candidate

arms <- list()
for(meth in c("quadratic", "conjugate")){
  for(tolv in tols){
    tag <- sprintf("pilot_%s_tol%s", if(meth == "quadratic") "quad" else "conjugate", tolv)
    message(sprintf("\n=== %s: method %s, tolPrefitLoops %.2f ===", tag, meth, tolv))
    res <- run_study(study = "study1", variant = tag, method = meth,
                     start_params = start_params, par_thomas = par_thomas,
                     rho_hat = rho_hat, K_hat = K_hat,
                     max.iter.prefit = 10, nRuns = nRuns, nCores = 32,
                     nSimsLevels = nSims, seedBases = seed,
                     tolPrefitLoops = tolv)
    arms[[tag]] <- list(method = meth, tol = tolv, fits = res[[1]]$res)
    save(arms, nRuns, nSims, seed, tols, start_params, truth, par_thomas, rho_hat, K_hat,
         file = paste0(resFolder, "study1_tolerance_pilot.RDa"))
  }
}

# Results: ----
L2         <- function(v) sqrt(sum((as.numeric(v) - truth)^2))
cycle_end  <- function(loop){ b <- Filter(Negate(is.null), loop); b[[length(b)]]$params }
finalErr   <- function(f) L2(f$params[nrow(f$params), ])
prefitErr  <- function(f) L2(cycle_end(f$prefitResults[[length(f$prefitResults)]]))
ok         <- function(fits) !vapply(fits, inherits, logical(1), "try-error")

cat("\n\nPer arm (medians over the fits that completed)\n")
cat("arm                        n   cycles  err@prefit   err@final   s/fit\n")
for(tag in names(arms)){
  a <- arms[[tag]]; f <- a$fits[ok(a$fits)]
  cat(sprintf("%-24s %3d    %4.1f     %.4f      %.4f    %5.0f\n", tag, length(f),
      median(vapply(f, function(x) length(x$prefitResults), numeric(1))),
      median(vapply(f, prefitErr, numeric(1))),
      median(vapply(f, finalErr,  numeric(1))),
      median(vapply(f, function(x) x$secs, numeric(1)))))
}

cat("\nPaired within method: 0.02 against 0.05, fit by fit\n")
for(meth in c("quadratic", "conjugate")){
  tagA <- sprintf("pilot_%s_tol0.05", if(meth == "quadratic") "quad" else "conjugate")
  tagB <- sprintf("pilot_%s_tol0.02", if(meth == "quadratic") "quad" else "conjugate")
  if(!all(c(tagA, tagB) %in% names(arms))) next
  A <- arms[[tagA]]$fits; B <- arms[[tagB]]$fits
  keep <- ok(A) & ok(B)                       # a pair is only usable if both arms ran
  eA <- vapply(A[keep], finalErr, numeric(1))
  eB <- vapply(B[keep], finalErr, numeric(1))
  tA <- vapply(A[keep], function(x) x$secs, numeric(1))
  tB <- vapply(B[keep], function(x) x$secs, numeric(1))
  w  <- suppressWarnings(wilcox.test(eA, eB, paired = TRUE))
  cat(sprintf("  %-10s %d pairs | median err %.4f -> %.4f (%+.1f%%) | 0.02 better in %d%% | p = %.3f | time %.0f -> %.0f s (%+.0f%%)\n",
      meth, sum(keep), median(eA), median(eB), 100*(median(eB)/median(eA) - 1),
      round(100*mean(eB < eA)), w$p.value, median(tA), median(tB),
      100*(median(tB)/median(tA) - 1)))
}
cat("\nA smaller tolerance buys cycles, and cycles cost time. Read the error change and\n")
cat("the time change together: this is a cost/accuracy trade, not a free improvement.\n")
