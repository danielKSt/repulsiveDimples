# Helpers for the trust region simulation studies.
#
# Sourced by results_analysis.R. Only ggplot2 needs to be attached; everything else is
# base R, so this can be sourced into a bare session without pulling in the tidyverse.

# Parameters at the end of each prefit loop and at the end of the main fit, for one fit.
#
# One row per stage: "1" .. "<n prefit loops>" for the state after each prefit loop
# (taken at the end of the loop's mu block, which runs last), then "final" for the main
# fit. Parameters come back on the natural scale.
#
# res   a single fit, as returned by min_contrast_trust_region
# nSims the ensemble size that fit used, carried through as a column
getParamsFinal <- function(res, nSims){
  finalPars <- exp(res$params[nrow(res$params), ])
  paramsRes <- data.frame(kappa = finalPars[1],
                          omega = finalPars[2],
                          mu    = finalPars[3],
                          prefitNumber = "final",
                          nSims = nSims,
                          row.names = NULL)

  # seq_along rather than 1:length, so a fit run with max.iter.prefit = 0 (prefitResults
  # is NULL) gives just the "final" row instead of erroring on 1:0.
  for(i in seq_along(res$prefitResults)){
    loopPars <- exp(res$prefitResults[[i]]$mu$params)
    # Named rather than positional, so this does not silently mis-assign if the columns
    # above are ever reordered.
    paramsRes <- rbind(paramsRes,
                       data.frame(kappa = loopPars[1],
                                  omega = loopPars[2],
                                  mu    = loopPars[3],
                                  prefitNumber = as.character(i),
                                  nSims = nSims,
                                  row.names = NULL))
  }
  return(paramsRes)
}

# The same thing over every fit at every ensemble size, stacked into one data frame.
#
# resList a list with one element per ensemble size, each of the form
#         list(res = <list of fits>, nSims = <ensemble size>), which is what
#         study<N>_run.R's run_level returns.
# quiet   set TRUE to suppress the note about fits that failed
#
# Fits that came back as "try-error" (mclapply hands those back rather than aborting the
# rest of the level) are dropped, and counted in a message so a silently thinned level
# does not pass unnoticed.
combineParamsFinal <- function(resList, quiet = FALSE){
  perLevel <- lapply(resList, function(level){
    failed <- vapply(level$res, inherits, logical(1), "try-error")
    if(!quiet && any(failed)){
      message(sprintf("nSims = %s: dropping %d of %d fits that returned an error",
                      level$nSims, sum(failed), length(failed)))
    }
    kept <- level$res[!failed]
    if(length(kept) == 0){
      return(NULL)
    }
    do.call(rbind, lapply(kept, getParamsFinal, nSims = level$nSims))
  })
  return(do.call(rbind, perLevel))
}

# Boxplots of the estimates, one panel per parameter and one box per (stage, nSims), so
# both the march towards the truth over the prefit loops and the narrowing spread with
# ensemble size are readable in one figure.
#
# paramsFinal   as returned by combineParamsFinal
# par_thomas    the true parameters, as stored in the study's data file
# start_params  starting parameters ON THE LOG SCALE, matching the run scripts; pass NULL
#               to leave the starting-value reference line off
# parameters    which parameters to show, and in what order
# log_scale     TRUE plots the y axes on log10, the scale the optimizer works on, where
#               equal relative spread reads as equal visual spread
plotParamsFinal <- function(paramsFinal, par_thomas, start_params = NULL,
                            parameters = c("kappa", "omega", "mu"),
                            log_scale = FALSE){
  # "final" has to be forced last: left as a character column it sorts ahead of the digits.
  # Levels come from the data so this survives a change to maxPrefitLoops, and a set of
  # fits that stopped pre-fitting on tolPrefitLoops after differing numbers of cycles.
  stageLevels <- c(sort(setdiff(unique(paramsFinal$prefitNumber), "final")), "final")

  long <- do.call(rbind, lapply(parameters, function(pn)
    data.frame(stage = factor(paramsFinal$prefitNumber, levels = stageLevels),
               nSims = factor(paramsFinal$nSims, levels = sort(unique(paramsFinal$nSims))),
               parameter = factor(pn, levels = parameters),
               value = paramsFinal[[pn]],
               row.names = NULL)))

  # One row per panel, so the reference lines survive faceting.
  refLines <- data.frame(parameter = factor(parameters, levels = parameters),
                         truth = as.numeric(unlist(par_thomas[parameters])))

  p <- ggplot2::ggplot(long, ggplot2::aes(x = stage, y = value, fill = nSims)) +
    ggplot2::geom_hline(data = refLines, ggplot2::aes(yintercept = truth),
                        colour = "red", linetype = 2)

  if(!is.null(start_params)){
    refLines$start <- as.numeric(exp(start_params)[match(parameters,
                                                         c("kappa", "omega", "mu"))])
    p <- p + ggplot2::geom_hline(data = refLines, ggplot2::aes(yintercept = start),
                                 colour = "purple", linetype = 3)
  }

  p <- p +
    ggplot2::geom_boxplot(outlier.size = 0.5, linewidth = 0.3,
                          position = ggplot2::position_dodge(width = 0.8)) +
    ggplot2::facet_wrap(~ parameter, ncol = 1, scales = "free_y") +
    ggplot2::scale_fill_brewer(palette = "Blues", name = "nSims") +
    ggplot2::labs(x = "prefit loop (then the main fit)", y = NULL,
                  subtitle = paste("red dashed = truth",
                                   if(!is.null(start_params)) ", purple dotted = starting value"
                                   else "")) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")

  if(log_scale){
    p <- p + ggplot2::scale_y_log10()
  }
  return(p)
}
