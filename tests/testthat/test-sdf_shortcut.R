# ==============================================================================
# File: test-sdf_shortcut.R
# Purpose: c_mu_shortcut() and c_sigma_shortcut(): within their documented error
#          bounds at sampled n and beyond the grid, the constants as printed, and
#          the kept copy of the former formulas. Blocks 10, 11, 13 and 14 of the
#          former test script R/se_imputed_test.R of the research repository
#          (commit db76e19), assertions and tolerances unchanged.

test_that("c_mu_shortcut() is within its documented error bounds at sampled n, both versions", {
  # Bounds = the errors of the constants AS PRINTED on the whole grid [5, 1001]
  # (the fit's output, both versions). The sampled n include the two ends of
  # the grid and n = 9, where the errors for the mean under S1 and S3 peak.
  bound <- list(blue     = c(S1 = 0.0034, S2 = 0.00046, S3 = 0.0033),
                wan_luo_shi = c(S1 = 0.0034, S2 = 0.00046, S3 = 0.0034))
  for (v in names(bound)) for (sc in names(bound[[v]])) for (nn in c(5L, 9L, 13L, 17L, 41L, 101L, 1001L)) {
    relerr <- abs(c_mu_shortcut(nn, sc, version = v) / c_mu(nn, sc, method = meth_mu[[v]]) - 1)
    expect_lte(relerr, bound[[v]][[sc]])
  }
})

test_that("c_sigma_shortcut() is within its documented error bounds at sampled n, both versions, with no upper limit on n", {
  # Bounds = the printed-constant errors of ?c_sigma_shortcut on the whole grid
  # [5, 1001] (same output file). c_sigma() consistency with
  # imputed_moments()$var_sd is checked in the same loop.
  bound_sd <- list(blue     = c(S1 = 0.0037, S2 = 0.0013, S3 = 0.0042),
                   wan_luo_shi = c(S1 = 0.0037, S2 = 0.00099, S3 = 0.0044))
  for (v in names(bound_sd)) for (sc in names(bound_sd[[v]])) for (nn in c(5L, 9L, 13L, 17L, 29L, 41L, 101L, 1001L)) {
    cs <- c_sigma(nn, sc, method = meth_sd(v, sc))
    expect_lt(abs(cs - sqrt(imputed_moments(nn, sc, NA, meth_sd(v, sc))$var_sd)), 1e-14)
    expect_lte(abs(c_sigma_shortcut(nn, sc, version = v) / cs - 1), bound_sd[[v]][[sc]])
  }
  # one formula serves every n >= 5, for c_mu_shortcut() as for c_sigma_shortcut(): no warning
  # at a large n, a warning below 5
  for (sc in c("S1", "S2", "S3")) for (v in c("blue", "wan_luo_shi")) {
    expect_false(warned(c_mu_shortcut(5000L, sc, version = v)))
    expect_false(warned(c_sigma_shortcut(5000L, sc, version = v)))
    expect_true(warned(c_mu_shortcut(3L, sc, version = v)))
    expect_true(warned(c_sigma_shortcut(3L, sc, version = v)))
  }
})

test_that("the shortcut formulas carry the constants as printed, both versions", {
  # Typed here from the companion paper's Table 3 (BLUE) and its supplement table of the
  # Wan-Luo-Shi version's formulas, so a reversion to the unrounded fit, or any slip
  # in a constant, fails this block (the bound tests above test accuracy bounds, which
  # a small slip would pass).
  nn <- c(5, 13, 41, 101, 1001); L <- log(nn)
  expect_lt(max(abs(c_mu_shortcut(nn, "S1") - 1.2533 / sqrt(nn + 2.006 * L - 0.810))), 1e-12)
  expect_lt(max(abs(c_mu_shortcut(nn, "S2") - 1.0780 / sqrt(nn + 0.142))), 1e-12)
  expect_lt(max(abs(c_mu_shortcut(nn, "S3") - 1.0780 / sqrt(nn + 0.838 * L - 0.576))), 1e-12)
  expect_lt(max(abs(c_sigma_shortcut(nn, "S1") - 0.8278 / (L + 0.330)^1.2041)), 1e-12)
  expect_lt(max(abs(c_sigma_shortcut(nn, "S2") - 1.1664 / sqrt(nn - 0.863))), 1e-12)
  expect_lt(max(abs(c_sigma_shortcut(nn, "S3") - 1.1664 / sqrt(nn + 5.267 * L^2 - 10.364 * L + 8.327))), 1e-12)
  expect_lt(max(abs(c_mu_shortcut(nn, "S1", "wan_luo_shi") - 1.2533 / sqrt(nn + 2.006 * L - 0.811))), 1e-12)
  expect_lt(max(abs(c_mu_shortcut(nn, "S2", "wan_luo_shi") - 1.0780 / sqrt(nn + 0.142))), 1e-12)
  expect_lt(max(abs(c_mu_shortcut(nn, "S3", "wan_luo_shi") - 1.0780 / sqrt(nn + 0.835 * L - 0.572))), 1e-12)
  expect_lt(max(abs(c_sigma_shortcut(nn, "S1", "wan_luo_shi") - 0.8705 / (L + 0.418)^1.2199)), 1e-12)
  expect_lt(max(abs(c_sigma_shortcut(nn, "S2", "wan_luo_shi") - 1.1664 / sqrt(nn - 0.845))), 1e-12)
  expect_lt(max(abs(c_sigma_shortcut(nn, "S3", "wan_luo_shi") - 1.1664 / sqrt(nn + 5.191 * L^2 - 10.204 * L + 8.525))), 1e-12)
})

test_that("the formulas stay within their documented errors at n = 2001, and the kept copy of the former formulas is as printed", {
  # (a) At n = 2001 the formulas stay within the errors ?c_mu_shortcut and
  #     ?c_sigma_shortcut document for the sample sizes beyond the grid (the bounds
  #     below are those figures at n = 2001, from 03-fits.csv, rounded up).
  bey <- list(blue     = list(mu = c(S1 = 0.0008, S2 = 0.00003, S3 = 0.0009), sd = c(S1 = 0.008, S2 = 0.00003, S3 = 0.0032)),
              wan_luo_shi = list(mu = c(S1 = 0.0006, S2 = 0.00003, S3 = 0.0008), sd = c(S1 = 0.008, S2 = 0.00002, S3 = 0.0027)))
  for (v in names(bey)) for (sc in c("S1", "S2", "S3")) {
    mo <- osmoments::os_moments(2001L, scenario_ranks(2001L, sc))
    expect_lte(abs(c_mu_shortcut(2001L, sc, version = v) / c_mu(2001L, sc, method = meth_mu[[v]], moments = mo) - 1), bey[[v]]$mu[[sc]])
    expect_lte(abs(c_sigma_shortcut(2001L, sc, version = v) / c_sigma(2001L, sc, method = meth_sd(v, sc), moments = mo) - 1), bey[[v]]$sd[[sc]])
  }
  # (b) .sdf_formula_former() carries the formulas as first printed (typed here
  #     from the earlier report and the former table), so
  #     that the studies run with them keep reproducing.
  f <- fivenumMeta:::.sdf_formula_former
  nn <- c(5, 13, 41, 101, 1001)
  expect_lt(max(abs(f(nn, "S1", "mu", "blue") - 0.9966 / (nn + 0.470)^0.4586)), 1e-12)
  expect_lt(max(abs(f(nn, "S3", "mu", "blue") - 0.9864 / (nn + 0.124)^0.4846)), 1e-12)
  expect_lt(max(abs(f(nn, "S1", "sd", "blue") - 0.9446 / (log(nn) + 0.478)^1.2662)), 1e-12)
  expect_lt(max(abs(f(nn, "S3", "sd", "blue") - 0.5667 / (nn - 2.073)^0.4052)), 1e-12)
  expect_lt(max(abs(f(nn, "S2", "mu", "wan_luo_shi") - 1.0780 / sqrt(nn + 0.142))), 1e-12)
  expect_lt(max(abs(f(nn, "S1", "mu", "wan_luo_shi") - 0.9973 / (nn + 0.475)^0.4588)), 1e-12)
  expect_lt(max(abs(f(nn, "S3", "mu", "wan_luo_shi") - 0.9821 / (nn + 0.076)^0.4837)), 1e-12)
  expect_lt(max(abs(f(nn, "S2", "sd", "wan_luo_shi") - 1.1664 / sqrt(nn - 0.845))), 1e-12)
  expect_lt(max(abs(f(nn, "S3", "sd", "wan_luo_shi") - 0.5671 / (nn - 1.971)^0.4047)), 1e-12)
  expect_lt(max(abs(f(nn, "S1", "sd", "wan_luo_shi") - 1.0040 / (log(nn) + 0.580)^1.2860)), 1e-12)
})
