# ==============================================================================
# File: test-sdf.R
# Purpose: imputed_moments(), c_mu(), c_sigma() and inflation_factors(): the
#          quadratic form against the entrywise expansion, the consistency of the
#          factors with imputed_moments() and the BLUE default, and the scenarios
#          S4 and S5. Blocks 8, 9 and 15 of the former test script
#          R/se_imputed_test.R of the research repository (commit db76e19),
#          assertions and tolerances unchanged.

test_that("the quadratic form w'Bw equals the entrywise (Luo-style) expansion", {
  # For an arbitrary asymmetric w, Var(w'X)/sigma^2 = w'Bw must equal the Luo-style
  # entrywise expansion sum_{r,s} w_r w_s E[Z_r Z_s] - (w'alpha)^2.
  mo25 <- osmoments::os_moments(25, scenario_ranks(25, "S3"))
  w <- c(0.10, 0.20, 0.30, 0.25, 0.15)
  im8 <- imputed_moments(25, "S3", w_mean = w, method_sd = NA, moments = mo25)
  M <- mo25$B + tcrossprod(mo25$alpha)
  expand <- sum(outer(w, w) * M) - sum(w * mo25$alpha)^2
  expect_lt(abs(im8$var_mean - expand), 1e-12)
})

test_that("c_mu() and c_sigma() are consistent with imputed_moments(), and their default method is the BLUE", {
  expect_lt(abs(c_mu(25, "S2") - sqrt(imputed_moments(25, "S2", "blue", NA)$var_mean)), 1e-14)
  expect_lt(abs(c_mu(25, "S2", "luo") - sqrt(imputed_moments(25, "S2", "luo", NA)$var_mean)), 1e-14)
  expect_lt(abs(c_sigma(25, "S2") - sqrt(imputed_moments(25, "S2", NA, "blue")$var_sd)), 1e-14)
})

test_that("S4 and S5: the ranks, the BLUE factors, inflation_factors() and se_imputed()", {
  # (a) the ranks at n = 21 (m = 1): the 10th and 90th percentiles at 2m+1 = 3 and
  #     18m+1 = 19, the quartiles at 6, 11, 16; S5 adds 1 and 21.
  # (b) the BLUE factors at n = 21 against the values of the companion paper's study
  #     of these scenarios (its Scenarios 5 and 7, written by the same engine), and inflation_factors()
  #     against them and against its definition C_mu sqrt(n), C_sigma / sqrt(1 - b_S^2).
  # (c) inflation_factors() under S2 at n = 41 against Table 2 of the paper (1.076,
  #     1.651), both versions consistent with c_mu() and c_sigma().
  # (d) se_imputed() under S4 applies the BLUE weights; the Wan-Luo-Shi version and
  #     presets, and a sample size below the scenario's smallest, are rejected.
  expect_identical(scenario_ranks(21L, "S4"), c(3L, 6L, 11L, 16L, 19L))
  expect_identical(scenario_ranks(21L, "S5"), c(1L, 3L, 6L, 11L, 16L, 19L, 21L))
  expect_identical(scenario_ranks(41L, "S4"), c(5L, 11L, 21L, 31L, 37L))
  pin <- rbind(S4 = c(c_mu = 0.224463076246977, c_sigma = 0.189701723919523,
                      kappa_mu = 1.02861903762442, kappa_sigma = 1.20753315338772),
               S5 = c(c_mu = 0.220903472809296, c_sigma = 0.162292532090282,
                      kappa_mu = 1.01230688542722, kappa_sigma = 1.03306184570781))
  bS21 <- sqrt(2 / 20) * exp(lgamma(21 / 2) - lgamma(20 / 2))
  for (sc in c("S4", "S5")) {
    cm <- c_mu(21L, sc); cs <- c_sigma(21L, sc); k <- inflation_factors(21L, sc)
    expect_lt(abs(cm - pin[sc, "c_mu"]), 1e-12)
    expect_lt(abs(cs - pin[sc, "c_sigma"]), 1e-12)
    expect_lt(abs(k[["kappa_mu"]] - pin[sc, "kappa_mu"]), 1e-12)
    expect_lt(abs(k[["kappa_sigma"]] - pin[sc, "kappa_sigma"]), 1e-12)
    expect_lt(abs(k[["kappa_mu"]] - cm * sqrt(21)), 1e-14)
    expect_lt(abs(k[["kappa_sigma"]] - cs / sqrt(1 - bS21^2)), 1e-12)
  }
  k2 <- inflation_factors(41L, "S2"); k2s <- inflation_factors(41L, "S2", version = "wan_luo_shi")
  expect_identical(round(unname(k2), 3), c(1.076, 1.651))
  expect_lt(abs(k2s[["kappa_mu"]] - c_mu(41L, "S2", method = "luo") * sqrt(41)), 1e-12)
  expect_lt(abs(k2s[["kappa_sigma"]] / k2[["kappa_sigma"]] -
                c_sigma(41L, "S2", method = "wan") / c_sigma(41L, "S2")), 1e-12)
  q7 <- c(min = 1, p10 = 6, q1 = 10, med = 12, q3 = 15, p90 = 20, max = 30)
  r4 <- se_imputed(q7, 21L, "S4"); mo4 <- osmoments::os_moments(21L, scenario_ranks(21L, "S4"))
  bw4 <- fivenumMeta:::.blue_weights(mo4$alpha, mo4$B); x4 <- as.numeric(q7[c("p10", "q1", "med", "q3", "p90")])
  expect_true(isTRUE(r4$exact)); expect_identical(r4$version_used, "blue")
  expect_lt(abs(r4$mu_hat - sum(bw4$w_mu * x4)), 1e-12)
  expect_lt(abs(r4$sigma_hat - sum(bw4$w_sig * x4)), 1e-12)
  expect_lt(abs(r4$se_mu - pin["S4", "c_mu"] * r4$sigma_hat), 1e-12)
  expect_true(bad(se_imputed(q7, 21L, "S4", version = "wan_luo_shi")))
  expect_true(bad(se_imputed(q7, 6L, "S5")))
  expect_true(bad(estimator_weights(21L, "S4", "mean", "luo")))
  expect_true(bad(inflation_factors(21L, "S5", version = "wan_luo_shi")))
  expect_true(bad(c_mu(4L, "S4")))
})

test_that(".b_s() and .var_s() agree with the gamma form at moderate n and stay in range at large n", {
  bS <- function(n) sqrt(2 / (n - 1)) * exp(lgamma(n / 2) - lgamma((n - 1) / 2))
  for (n in c(2, 3, 21, 100, 1001)) {
    expect_lt(abs(fivenumMeta:::.b_s(n) - bS(n)), 1e-12)
    expect_lt(abs(fivenumMeta:::.var_s(n) - (1 - bS(n)^2)), 1e-12)
  }
  expect_lt(abs(fivenumMeta:::.b_s(2) - sqrt(2 / pi)), 1e-15)
  # the lgamma difference gave b_S > 1 at n = 1e8; the beta form stays below 1 and Var(S) near 1/(2n)
  for (n in c(1e6, 1e8, 1e9)) {
    expect_lt(fivenumMeta:::.b_s(n), 1)
    expect_lt(abs(fivenumMeta:::.var_s(n) / (1 / (2 * n)) - 1), 1e-5)
  }
  expect_true(is.finite(inflation_factors(1001L, "S2")[["kappa_sigma"]]))
  # n outside the integer range, or infinite, is refused with a message, not a coercion warning
  expect_error(c_mu(2^31, "S2"), "integer")
  expect_error(c_mu(Inf, "S2"), "integer")
  expect_error(estimator_weights(2^31, "S2", "mean", "luo"), "integer")
  expect_error(scenario_ranks(2^31, "S2"), "integer")
})

test_that("imputed_moments() defaults to the BLUE of both, and method_sd = NULL keeps the Wan-Luo-Shi pairing", {
  im <- imputed_moments(25, "S3")
  expect_equal(im$w_mean, estimator_weights(25, "S3", "mean", "blue"), tolerance = 1e-15)
  expect_equal(im$w_sd, estimator_weights(25, "S3", "sd", "blue"), tolerance = 1e-15)
  expect_lt(abs(im$c_n - c_mu(25, "S3")), 1e-15)
  expect_lt(abs(sqrt(im$var_sd) - c_sigma(25, "S3")), 1e-15)
  # the Wan-Luo-Shi version: Luo's mean with method_sd = NULL, Shi's SD under S3 and Wan's under S2
  expect_equal(imputed_moments(25, "S3", "luo", NULL)$w_sd, estimator_weights(25, "S3", "sd", "shi"), tolerance = 1e-15)
  expect_equal(imputed_moments(25, "S2", "luo", NULL)$w_sd, estimator_weights(25, "S2", "sd", "wan"), tolerance = 1e-15)
  expect_equal(imputed_moments(25, "S2", "luo", NULL)$w_mean, estimator_weights(25, "S2", "mean", "luo"), tolerance = 1e-15)
  expect_lt(abs(inflation_factors(25, "S3", version = "wan_luo_shi")[["kappa_mu"]] -
                imputed_moments(25, "S3", "luo", NULL)$c_n * sqrt(25)), 1e-14)
  # NA skips a part, as before
  expect_true(is.null(imputed_moments(25, "S2", NA, "blue")$w_mean))
})
