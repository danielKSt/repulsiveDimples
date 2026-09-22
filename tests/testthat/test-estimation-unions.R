# The union estimators place each independent snapshot in its own box of a disjoint
# union window, so that one spatstat call over the union estimates the same summary
# function the snapshots share. Both loops that build that union used to be written
# with `:`, which counts backwards when the sequence is empty:
#
#   1:(length(snapshots)-1)        becomes 1:0 == c(1, 0) for a single snapshot,
#                                  unioning in a spurious empty box, and
#   (t_first+1):length(snapshots)  becomes 2:1 == c(2, 1), which reads snapshots[2]
#                                  (NA, silently NULL) and then re-adds snapshot 1
#                                  at offset zero, duplicating every point.
#
# Nothing errored: `points_input[[NA]]` is NULL, `NULL[, 1:2] + offset` is numeric(0),
# and `rbind(df, numeric(0))` is a no-op. Only the numbers were wrong -- on study1's
# parameters K(1.0) came back as 9.148 instead of 3.824. A single snapshot is a
# legitimate real-data case, so these tests pin it down.

sidelength <- 40
r_vec <- seq(from = 0, to = 3, by = 0.05)

# study1's parameters, the configuration the bug was measured on.
one_pattern <- function(seed = 99) {
  set.seed(seed)
  sim <- rThomas_matern_thinned(kappa = 0.1234, scale = 0.8, mu = 16.2075,
                                repulsionRange = 0.4,
                                xlims = c(0, sidelength), ylims = c(0, sidelength),
                                saveparents = TRUE)
  sim$thinned
}

as_ppp <- function(pat, l = sidelength) {
  spatstat.geom::ppp(x = pat$x, y = pat$y,
                     window = spatstat.geom::owin(c(0, l), c(0, l)))
}


test_that("a single snapshot reproduces spatstat's Kest on that pattern", {
  pat <- one_pattern()

  union_K <- K_est.unions(points_input = list(pat), l = sidelength, spacing = 5,
                          r_vec = r_vec, timescale = 1)
  direct_K <- spatstat.explore::Kest(as_ppp(pat), r = r_vec, correction = "border")

  # One snapshot means one box holding one copy of the pattern, so this is not merely
  # close to the direct estimate -- it is the same computation.
  expect_equal(union_K$border, direct_K$border, tolerance = 1e-12)
  expect_equal(union_K$r, direct_K$r, tolerance = 1e-12)
})

test_that("a single snapshot is not silently duplicated", {
  pat <- one_pattern()

  # The duplicate copy landed on top of the original, so spatstat warned about
  # duplicated points and roughly doubled every pair count.
  expect_no_warning(
    K_est.unions(points_input = list(pat), l = sidelength, spacing = 5,
                 r_vec = r_vec, timescale = 1)
  )
})

test_that("a single snapshot does not union in an empty second box", {
  pat <- one_pattern()

  # The spurious box doubled the window area without adding points. Border-corrected
  # K divides by the intensity estimated over the whole union, so an extra empty box
  # inflates K even once the duplication is gone. Comparing against Kest on the plain
  # window pins the area down; comparing against the pre-fix value guards the number
  # that was actually observed.
  union_K <- K_est.unions(points_input = list(pat), l = sidelength, spacing = 5,
                          r_vec = r_vec, timescale = 1)
  at_one <- union_K$border[which.min(abs(union_K$r - 1))]

  expect_equal(at_one,
               spatstat.explore::Kest(as_ppp(pat), r = r_vec,
                                      correction = "border")$border[which.min(abs(r_vec - 1))],
               tolerance = 1e-12)
  expect_lt(at_one, 5)   # the broken version returned 9.148 here
})

test_that("m identical snapshots give the same K as a single snapshot", {
  pat <- one_pattern()

  one <- K_est.unions(points_input = list(pat), l = sidelength, spacing = 5,
                      r_vec = r_vec, timescale = 1)

  for (m in c(2, 5)) {
    many <- K_est.unions(points_input = rep(list(pat), m), l = sidelength, spacing = 5,
                         r_vec = r_vec, timescale = 1)
    # The boxes are disjoint and spacing*l apart, so no pair crosses between copies and
    # each copy's border distances are its own box's. The estimate is therefore the
    # same up to the arithmetic of the larger window.
    expect_equal(many$border, one$border, tolerance = 1e-6)
  }
})

test_that("K_est.unions handles snapshots whose only pattern is the last one", {
  pat <- one_pattern()

  # Same backwards loop, reached a different way: leading empty snapshots push t_first
  # up to length(snapshots), so (t_first+1):length(snapshots) counted backwards and
  # re-added the pattern displaced into another box. Both orderings describe one
  # pattern observed in one of three windows, so they must agree.
  trailing <- K_est.unions(points_input = list(0, 0, pat), l = sidelength, spacing = 5,
                           r_vec = r_vec, timescale = 1)
  leading  <- K_est.unions(points_input = list(pat, 0, 0), l = sidelength, spacing = 5,
                           r_vec = r_vec, timescale = 1)

  expect_equal(trailing$border, leading$border, tolerance = 1e-12)
})

test_that("pcf_est.unions takes its points in the same coordinates as K_est.unions", {
  # pcf_est.unions used to scale its input by l while K_est.unions did not, so the two
  # estimators disagreed about what a point pattern was. An ordinary observation window
  # pattern handed to the old version landed at l times its coordinates, far outside the
  # window, and spatstat discarded all of it.
  pat <- one_pattern()
  r_grid <- seq(from = 0, to = 3, by = 0.05)

  union_pcf <- suppressWarnings(
    pcf_est.unions(points_input = list(pat), bw = 0.15, l = sidelength,
                   spacing = 5, rMax = 3, dr = 0.05, timescale = 1))
  direct_pcf <- spatstat.explore::pcf.ppp(as_ppp(pat), r = r_grid, bw = 0.15,
                                          divisor = "r", zerocor = NULL)

  # pcf.ppp returns its estimate in `iso`; guard against comparing nothing at all.
  expect_length(union_pcf$iso, length(r_grid))
  expect_true(any(is.finite(union_pcf$iso)))

  expect_equal(union_pcf$iso, direct_pcf$iso, tolerance = 1e-12)
})

test_that("pcf_est.unions keeps its points inside the observation window", {
  # Separate block: a failing expect_no_warning aborts the rest of its test_that, which
  # would mask the comparison above.
  pat <- one_pattern()

  expect_no_warning(
    pcf_est.unions(points_input = list(pat), bw = 0.15, l = sidelength,
                   spacing = 5, rMax = 3, dr = 0.05, timescale = 1)
  )
})

test_that("both estimators lay snapshots out on the same union window", {
  # Rebuild the union by hand and check each estimator against spatstat applied to it.
  # This pins the shared geometry -- box positions, per-snapshot offsets and the
  # coordinate convention -- rather than only the single-snapshot special case.
  pat <- one_pattern()
  m <- 3; spacing <- 5

  win <- spatstat.geom::owin(c(0, sidelength), c(0, sidelength))
  for (t_ind in seq_len(m - 1)) {
    win <- spatstat.geom::union.owin(
      win,
      spatstat.geom::owin(xrange = c(spacing*t_ind*sidelength, (spacing*t_ind+1)*sidelength),
                          yrange = c(spacing*t_ind*sidelength, (spacing*t_ind+1)*sidelength)))
  }
  pts <- pat[, 1:2]
  for (t_ind in 2:m) pts <- rbind(pts, pat[, 1:2] + spacing*(t_ind-1)*sidelength)
  by_hand <- spatstat.geom::ppp(x = pts$x, y = pts$y, window = win)

  inp <- rep(list(pat), m)

  expect_equal(
    K_est.unions(points_input = inp, l = sidelength, spacing = spacing,
                 r_vec = r_vec, timescale = 1)$border,
    spatstat.explore::Kest(by_hand, r = r_vec, correction = "border")$border,
    tolerance = 1e-12)

  got_pcf <- suppressWarnings(
    pcf_est.unions(points_input = inp, bw = 0.15, l = sidelength, spacing = spacing,
                   rMax = 3, dr = 0.05, timescale = 1))
  expect_true(any(is.finite(got_pcf$iso)))
  expect_equal(
    got_pcf$iso,
    spatstat.explore::pcf.ppp(by_hand, r = seq(from = 0, to = 3, by = 0.05),
                              bw = 0.15, divisor = "r", zerocor = NULL)$iso,
    tolerance = 1e-12)
})

test_that("m identical snapshots leave the pcf unchanged in shape", {
  pat <- one_pattern()
  n <- nrow(pat)
  f <- function(m) pcf_est.unions(points_input = rep(list(pat), m), bw = 0.15,
                                  l = sidelength, spacing = 5, rMax = 3, dr = 0.05,
                                  timescale = 1)
  one <- f(1)

  for (m in c(2, 3, 5)) {
    many <- f(m)
    ok <- is.finite(one$iso) & is.finite(many$iso) & one$iso > 1e-8
    expect_true(sum(ok) > 50)

    # Unlike Kest, pcf.ppp normalises by the unbiased lambda^2 = n(n-1)/|W|^2. Stacking
    # m copies leaves n(n-1) growing slightly faster than m^2, so the estimate is
    # rescaled by exactly m(n-1)/(mn-1) -- about 4e-4 here. That is a constant in r: the
    # shape of the estimate is untouched, which is the property worth pinning down.
    ratio <- many$iso[ok] / one$iso[ok]
    expect_equal(max(ratio) - min(ratio), 0, tolerance = 1e-7)
    expect_equal(mean(ratio), m*(n - 1)/(m*n - 1), tolerance = 1e-6)

    # and the values themselves still agree to that same small factor
    expect_equal(many$iso[ok], one$iso[ok], tolerance = 1e-3)
  }
})

test_that("intensity_est is unaffected by the number of snapshots", {
  pat <- one_pattern()

  # intensity_est loops over 1:length(snapshots), which is well behaved for one
  # snapshot; this pins that down alongside the estimators that were not.
  one <- intensity_est(points_input = list(pat), l = sidelength, timescale = 1)
  expect_equal(one, nrow(pat) / sidelength^2, tolerance = 1e-12)

  for (m in c(2, 5)) {
    expect_equal(intensity_est(points_input = rep(list(pat), m), l = sidelength,
                               timescale = 1),
                 one, tolerance = 1e-12)
  }
})
