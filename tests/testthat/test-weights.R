# ==============================================================================
# File: test-weights.R
# Purpose: The BLUE's structure (Balakrishnan et al. 2022) and the preset weights'
#          identities. Blocks 6 and 7 of the former test script R/se_imputed_test.R
#          of the research repository (commit db76e19), assertions and tolerances
#          unchanged, with the moments from osmoments::os_moments().

test_that("the BLUE has the structure of Balakrishnan et al. (2022): scale = range / E[range] under S1 and S2, zero covariance, unbiased", {
  # S1/S2 scale BLUE is Range / E[Range]: weights (-1, 0, 1) / (2 alpha_max)
  # (eqs. 30-32); the symmetric-case covariance is zero (eq. 97).
  for (sc in c("S1", "S2")) {
    mo <- osmoments::os_moments(41, scenario_ranks(41, sc))
    wsd <- estimator_weights(41, sc, "sd", "blue", moments = mo)
    expect_lt(wsd %~% (c(-1, 0, 1) / (2 * mo$alpha[3])), 1e-7)
  }
  im <- imputed_moments(41, "S3", method_mean = "blue", method_sd = "blue")
  expect_lt(abs(im$cov), 1e-8)
  # BLUE unbiasedness: w_mu'1 = 1, w_mu'alpha = 0; w_sd'1 = 0, w_sd'alpha = 1.
  expect_lt(abs(im$gaps) %~% 0, 1e-7)
})

test_that("the preset weights satisfy their sum constraints, and the Blom gap is bounded", {
  # Symmetric mean presets: w'1 = 1 exactly and w'alpha = 0 to integration error.
  # Wan/Shi SD presets: w'1 = 0 exactly; w'alpha - 1 is the (small) Blom gap.
  mo25 <- osmoments::os_moments(25, scenario_ranks(25, "S3"))
  for (meth in c("wan", "luo")) {
    w <- estimator_weights(25, "S3", "mean", meth)
    expect_lt(abs(sum(w) - 1), 1e-12)
    expect_lt(abs(sum(w * mo25$alpha)), 1e-9)
  }
  for (meth in c("wan", "shi")) {
    w <- estimator_weights(25, "S3", "sd", meth)
    expect_lt(abs(sum(w)), 1e-12)
    expect_lt(abs(sum(w * mo25$alpha) - 1), 0.02)
  }
  # Luo/Shi weights at n = 25 agree with hand-checked values (same formulas).
  w1 <- 2.2 / (2.2 + 25^0.75); w2 <- 0.7 - 0.72 / 25^0.55
  expect_lt(estimator_weights(25, "S3", "mean", "luo") %~% c(w1 / 2, w2 / 2, 1 - w1 - w2, w2 / 2, w1 / 2), 1e-15)
})

test_that("a supplied moments object is checked for its ranks, n, distribution and shape; a bare list for its size", {
  r <- scenario_ranks(25, "S2"); good <- osmoments::os_moments(25, r); ref <- c_mu(25, "S2")
  expect_identical(c_mu(25, "S2", moments = good), ref)
  expect_lt(abs(c_mu(25, "S2", moments = good[c("alpha", "B")]) - ref), 1e-15)
  expect_error(c_mu(25, "S2", moments = osmoments::os_moments(30, r)), "n = 30")
  expect_error(c_mu(25, "S2", moments = osmoments::os_moments(25, r, dist = "unif")), "normal")
  expect_error(c_mu(25, "S2", moments = osmoments::os_moments(25, c(1, 13, 25))), "ranks")
  expect_error(c_mu(25, "S2", moments = osmoments::os_moments(25, r, location = 5)), "location 0")
  bad <- good; bad$B[3, 1] <- 10 * bad$B[3, 1]
  expect_error(c_mu(25, "S2", moments = bad), "symmetric")
  bad <- good; bad$alpha[2] <- NA
  expect_error(c_mu(25, "S2", moments = bad), "finite")
  expect_error(c_mu(25, "S2", moments = list(B = good$B)), "alpha and B")
})

test_that("a singular design is refused for the BLUE alone, with the package's message, and fixed weights still have a variance", {
  # under type 3 at n = 3 the third quartile is the median: no BLUE, but Luo's weights have a variance
  expect_error(c_mu(3, "S2", method = "blue", type = 3), "linearly dependent")
  expect_error(estimator_weights(3, "S2", "mean", "blue", type = 3), "linearly dependent")
  des <- fivenumMeta:::.reported_design(3, "S2", 3); mo <- osmoments::os_moments(3, c(1, 2))
  w <- estimator_weights(3, "S2", "mean", "luo")
  expect_lt(abs(c_mu(3, "S2", method = "luo", type = 3) - sqrt(drop(w %*% t(des$W) %*% mo$B %*% des$W %*% w))), 1e-14)
  expect_false(fivenumMeta:::.reported_moments(3, "S2", 3)$full_rank)
  expect_true(fivenumMeta:::.reported_moments(3, "S2", 7)$full_rank)
  # the four designs where the Cholesky test let a singular matrix through to solve()
  for (case in list(c(10L, 1L), c(10L, 3L), c(10L, 4L), c(13L, 3L)))
    expect_error(c_mu(case[1], "S5", type = case[2]), "linearly dependent")
  expect_error(fivenumMeta:::.blue_weights(c(-1, 0, 1), matrix(1, 3, 3)), "singular")
})

test_that("the engine's moments are kept for the session, and the quantile columns of W sum to one", {
  r <- scenario_ranks(25, "S2"); key <- paste(25, paste(r, collapse = ","))
  rm(list = intersect(key, ls(fivenumMeta:::.moments_cache)), envir = fivenumMeta:::.moments_cache)
  m1 <- fivenumMeta:::.reported_moments(25, "S2", 7)$moments
  expect_true(exists(key, envir = fivenumMeta:::.moments_cache, inherits = FALSE))
  m2 <- fivenumMeta:::.reported_moments(25, "S2", 7)$moments
  expect_identical(m1, m2)
  expect_equal(m1$alpha, osmoments::os_moments(25, r)$alpha)
  # a residue dropped under type 8 at n = 41 leaves the column summing to one exactly
  des <- fivenumMeta:::.reported_design(41L, "S2", 8)
  expect_lt(max(abs(colSums(des$W) - 1)), 1e-15)
})

test_that("the default method of estimator_weights() is the BLUE, for both estimands and every scenario", {
  # gives the BLUE, so it works for the SD too (the former default "luo" had no SD preset)
  for (sc in c("S1", "S2", "S3")) for (est in c("mean", "sd"))
    expect_identical(estimator_weights(25, sc, est), estimator_weights(25, sc, est, "blue"))
  expect_identical(estimator_weights(21, "S4", "sd"), estimator_weights(21, "S4", "sd", "blue"))
  expect_identical(estimator_weights(21, "S5", "mean"), estimator_weights(21, "S5", "mean", "blue"))
  expect_identical(estimator_weights(30, "S2", "mean"), estimator_weights(30, "S2", "mean", "blue"))
})
