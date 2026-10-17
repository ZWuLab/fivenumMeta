# ==============================================================================
# File: test-se_imputed.R
# Purpose: se_imputed(): both versions, both standard errors, exact at every n, and
#          the input checks. Block 12 of the former test script R/se_imputed_test.R
#          of the research repository (commit db76e19), assertions and tolerances
#          unchanged, with the moments from osmoments::os_moments().

test_that("se_imputed() matches the published formulas (wan_luo_shi) and the BLUE weights (default), both SEs, exact at every n, and the input checks", {
  # (a) version = "wan_luo_shi": mu-hat, sigma-hat, the naive SE and both plug-in SEs
  #     against the source equations typed here independently of estimator_weights():
  #     Luo 2018 eqs. (7), (11), (15) for the mean; Wan 2014 eqs. (9), (16) and
  #     Shi 2020 eq. (9) for the SD, with Blom divisors. Checked at n = 25 (exact
  #     rank) and n = 30 (interpolated quartiles, type 7): the factors are c_mu() and
  #     c_sigma() of the same estimators at that n.
  # (b) version = "blue" (the default): the weights are Lloyd's from the means and
  #     covariance matrix of the reported values (checked against .blue_weights()
  #     applied to os_moments() at the scenario's ranks at n = 25, and to the moments
  #     of the interpolated values at n = 30), the factors are c_mu()/c_sigma() with
  #     method "blue", and the point estimates are those weights applied to q.
  q5 <- c(min = 3, q1 = 10, med = 12, q3 = 15, max = 30)
  a <- 3; q1 <- 10; m <- 12; q3 <- 15; b <- 30
  hand <- function(n, sc) {
    xi  <- 2 * qnorm((n - 0.375) / (n + 0.25))          # Wan eq. 9 divisor
    eta <- 2 * qnorm((0.75 * n - 0.125) / (n + 0.25))   # Wan eq. 16 divisor
    switch(sc,
      S1 = { w1 <- 4 / (4 + n^0.75)                                       # Luo eq. 7
             c(mu = w1 * (a + b) / 2 + (1 - w1) * m, sd = (b - a) / xi) },
      S2 = { w1 <- 0.7 + 0.39 / n                                         # Luo eq. 11
             c(mu = w1 * (q1 + q3) / 2 + (1 - w1) * m, sd = (q3 - q1) / eta) },
      S3 = { w1 <- 2.2 / (2.2 + n^0.75); w2 <- 0.7 - 0.72 / n^0.55       # Luo eq. 15
             th <- 1 / (1 + 0.07 * n^0.6)                                  # Shi eq. 9
             c(mu = w1 * (a + b) / 2 + w2 * (q1 + q3) / 2 + (1 - w1 - w2) * m,
               sd = th * (b - a) / xi + (1 - th) * (q3 - q1) / eta) })
  }
  xq <- function(sc) as.numeric(q5[switch(sc, S1 = c("min", "med", "max"), S2 = c("q1", "med", "q3"), S3 = names(q5))])
  for (sc in c("S1", "S2", "S3")) for (nn in c(25L, 30L)) {
    # (a) the Wan-Luo-Shi version, against the hand-typed source equations
    r <- se_imputed(q5, nn, sc, version = "wan_luo_shi"); h <- hand(nn, sc)
    cn <- c_mu(nn, sc, method = "luo")
    cs <- c_sigma(nn, sc, method = if (sc == "S3") "shi" else "wan")
    expect_true(isTRUE(r$exact)); expect_identical(r$version_used, "wan_luo_shi")
    expect_lt(abs(r$mu_hat - h[["mu"]]), 1e-12)
    expect_lt(abs(r$sigma_hat - h[["sd"]]), 1e-12)
    expect_lt(abs(r$se_naive - h[["sd"]] / sqrt(nn)), 1e-12)
    expect_lt(abs(r$se_mu - cn * h[["sd"]]), 1e-12)
    expect_lt(abs(r$se_sigma - cs * h[["sd"]]), 1e-12)
    expect_lt(abs(r$c_n - cn), 1e-12); expect_lt(abs(r$c_sigma - cs), 1e-12)
    # (b) the BLUE, the default
    rb <- se_imputed(q5, nn, sc)
    rm <- if (nn == 25L) { mo <- osmoments::os_moments(nn, scenario_ranks(nn, sc)); list(alpha = mo$alpha, B = mo$B) } else
      fivenumMeta:::.reported_moments(nn, sc, 7)
    bw <- fivenumMeta:::.blue_weights(rm$alpha, rm$B)
    expect_identical(rb$version_used, "blue"); expect_true(isTRUE(rb$exact))
    expect_lt(max(abs(rb$w_mean - bw$w_mu)), 1e-12); expect_lt(max(abs(rb$w_sd - bw$w_sig)), 1e-12)
    expect_lt(abs(rb$mu_hat - sum(bw$w_mu * xq(sc))), 1e-12)
    expect_lt(abs(rb$sigma_hat - sum(bw$w_sig * xq(sc))), 1e-12)
    expect_lt(abs(rb$c_n - c_mu(nn, sc, method = "blue")), 1e-12)
    expect_lt(abs(rb$c_sigma - c_sigma(nn, sc, method = "blue")), 1e-12)
    expect_lt(abs(rb$se_mu - rb$c_n * rb$sigma_hat), 1e-12)
  }
  # (c) S1 at an odd n not of the form 4Q+1 reports single order statistics: the BLUE
  #     equals the one at the ranks of scenario_ranks()
  mo7 <- osmoments::os_moments(7L, scenario_ranks(7L, "S1"))
  r7 <- se_imputed(q5, 7L, "S1"); bw7 <- fivenumMeta:::.blue_weights(mo7$alpha, mo7$B)
  expect_true(isTRUE(r7$exact)); expect_lt(max(abs(r7$w_mean - bw7$w_mu)), 1e-12)
  # (d) invalid inputs are rejected: unordered values, a missing entry, a factor
  #     (its level codes must not be taken as values), a wrong scenario, a
  #     non-integer n, an infinite n
  qf <- factor(c(10, 12, 15)); names(qf) <- c("q1", "med", "q3")
  expect_true(bad(se_imputed(c(q1 = 15, med = 12, q3 = 10), 25L, "S2")))
  expect_true(bad(se_imputed(c(q1 = 1, med = 2), 25L, "S2")))
  expect_true(bad(se_imputed(qf, 25L, "S2")))
  expect_true(bad(se_imputed(c(q1 = 10, med = 12, q3 = 15), 25L, "S6")))
  expect_true(bad(se_imputed(c(q1 = 10, med = 12, q3 = 15), 25.5, "S2")))
  expect_true(bad(se_imputed(c(q1 = 10, med = 12, q3 = 15), Inf, "S2")))
  expect_true(bad(se_imputed(c(q1 = 10, med = 12, q3 = 15), 25L, "S2", type = 10)))
})

test_that("a repeated required name in q, and n above the integer range, are refused", {
  expect_error(se_imputed(c(q1 = 10, q1 = 99, med = 12, q3 = 15), 25, "S2"), "more than once")
  expect_error(se_imputed(c(q1 = 10, med = 12, q3 = 15, med = 13), 25, "S2"), "med")
  # an extra, unrelated name is still ignored
  expect_lt(abs(se_imputed(c(q1 = 10, med = 12, q3 = 15, note = 1), 25, "S2")$mu_hat -
                se_imputed(c(q1 = 10, med = 12, q3 = 15), 25, "S2")$mu_hat), 1e-15)
  expect_error(se_imputed(c(q1 = 10, med = 12, q3 = 15), 2^31, "S2"), "integer")
  expect_false(warned(try(se_imputed(c(q1 = 10, med = 12, q3 = 15), 2^31, "S2"), silent = TRUE)))
})

test_that("the result has the class se_imputed with a compact print, and a constant summary gives an SD of exactly zero with a warning", {
  r <- se_imputed(c(q1 = 10, med = 12, q3 = 15), 25, "S2")
  expect_s3_class(r, "se_imputed"); expect_true(is.list(r))
  expect_identical(class(r[c("mu_hat", "se_mu")]), "list")       # subsetting gives a plain list
  expect_identical(nrow(as.data.frame(r[c("mu_hat", "sigma_hat", "se_mu")])), 1L)
  expect_output(print(r), "Imputed mean and SD under S2, n = 25, the BLUE, quantile type 7", fixed = TRUE)
  expect_output(print(r), "mu_hat"); expect_output(print(r), "C_sigma(n) =", fixed = TRUE)
  expect_output(print(se_imputed(c(q1 = 10, med = 12, q3 = 15), 25, "S2", version = "wan_luo_shi")),
                "the Wan-Luo-Shi estimators", fixed = TRUE)
  # a constant summary: the SD weights sum to zero, so sigma_hat is exactly zero, with a warning
  expect_warning(rc <- se_imputed(c(q1 = 5, med = 5, q3 = 5), 25, "S2"), "all equal")
  expect_identical(rc$sigma_hat, 0); expect_identical(rc$se_naive, 0)
  expect_identical(rc$se_mu, 0); expect_identical(rc$se_sigma, 0)
  expect_lt(abs(rc$mu_hat - 5), 1e-12)
  expect_warning(se_imputed(c(min = 2, q1 = 2, med = 2, q3 = 2, max = 2), 30, "S3", version = "wan_luo_shi"), "all equal")
})
