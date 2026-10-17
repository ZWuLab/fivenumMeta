# ==============================================================================
# File: test-exact-any-n.R
#          section 4.3): (1) at the exact-rank sizes under type 7 the reported values
#          are the order statistics of scenario_ranks(), so the new path reproduces
#          the former calculation there; (2) at even n under S1 the BLUE equals an
#          independent calculation (the median as the average of the two middle
#          order statistics); (3) at an n not of the form 4Q + 1 a simulation agrees
#          with the exact factors; (4) the quantile definition is checked, and a
#          rounding residue of the quantile weights adds no rank.

# the former calculation: the BLUE and its factors from os_moments() at the scenario's ranks
old_path <- function(n, sc) {
  mo <- osmoments::os_moments(n, scenario_ranks(n, sc))
  bw <- fivenumMeta:::.blue_weights(mo$alpha, mo$B)
  c(c_mu = sqrt(drop(bw$w_mu %*% mo$B %*% bw$w_mu)),
    c_sigma = sqrt(drop(bw$w_sig %*% mo$B %*% bw$w_sig)))
}

test_that("at n = 4Q + 1 (and 20m + 1) under type 7 the factors are those of the former calculation", {
  skip_on_cran()                                         # 39 sample sizes x 3 scenarios
  # the 39-point grid of the companion paper: every n = 4Q + 1 from 5 to 101 and 14 larger ones
  grid <- c(seq(5L, 101L, by = 4L), 121L, 141L, 161L, 181L, 201L, 253L, 301L, 353L, 401L,
            501L, 601L, 701L, 801L, 1001L)
  expect_length(grid, 39L)
  for (sc in c("S1", "S2", "S3")) for (n in grid) {
    o <- old_path(n, sc)
    mo <- fivenumMeta:::.reported_moments(n, sc, 7)
    expect_identical(mo$ranks, scenario_ranks(n, sc))
    expect_lt(abs(c_mu(n, sc) - o[["c_mu"]]), 1e-10)
    expect_lt(abs(c_sigma(n, sc) - o[["c_sigma"]]), 1e-10)
  }
  for (sc in c("S4", "S5")) for (n in c(21L, 41L, 101L)) {
    o <- old_path(n, sc)
    expect_lt(abs(c_mu(n, sc) - o[["c_mu"]]), 1e-10)
    expect_lt(abs(c_sigma(n, sc) - o[["c_sigma"]]), 1e-10)
  }
})

test_that("under S1 at n = 40, 38 and 15 the BLUE and its SEs equal an independent calculation", {
  # reference values from a separate calculation written for the companion paper: the BLUE
  # from (min, median, max) by generalized least squares on the three statistics, the median
  # at even n the average of the two middle order statistics, computed with osmoments
  # directly and not through this package; the data are the vitamin D studies of Nnoaham
  # and Clarke (2008)
  ref <- data.frame(
    n   = c(40, 40, 40, 38, 15, 15),
    min = c(2.25, 9.0, 43.75, 48.5, 16.75, 26.25),
    med = c(16.0, 27.25, 65.75, 69.5, 39.75, 65.5),
    max = c(74.25, 132.5, 130.5, 125, 89.25, 114.75),
    mu  = c(20.3710977970914, 35.7957417606057, 69.9492006927114, 72.9890682355836, 44.3069729993174, 67.2196124525726),
    se_mu = c(3.02513374001294, 5.18894467904998, 3.64486599925171, 3.32043513379091, 5.9185921995893, 7.22476427122279),
    sd  = c(16.6606720783702, 28.577680578871, 20.0737958722029, 17.8730679506387, 20.8823775780084, 25.4909022848792),
    se_sd = c(2.57988389739465, 4.42521751844777, 3.10840177915258, 2.81132847296772, 4.5484677391046, 5.55226751601044))
  rel <- function(x, y) abs(x / y - 1)
  for (i in seq_len(nrow(ref))) {
    r <- se_imputed(c(min = ref$min[i], med = ref$med[i], max = ref$max[i]), ref$n[i], "S1")
    expect_lt(rel(r$mu_hat, ref$mu[i]), 1e-10)
    expect_lt(rel(r$sigma_hat, ref$sd[i]), 1e-10)
    expect_lt(rel(r$se_mu, ref$se_mu[i]), 1e-10)
    expect_lt(rel(r$se_sigma, ref$se_sd[i]), 1e-10)
  }
})

test_that("a simulation at n = 40 under S2 (interpolated quartiles, type 7) agrees with the exact factors", {
  skip_on_cran()
  set.seed(20261015)
  R <- 1e5L; n <- 40L
  x <- matrix(rnorm(R * n), R, n)
  x <- t(apply(x, 1, sort))
  # the type-7 quartiles of each sample, from the weights of the reported values
  des <- fivenumMeta:::.reported_design(n, "S2", 7)
  Y <- x[, des$ranks] %*% des$W
  w_mu <- estimator_weights(n, "S2", "mean", "blue"); w_sd <- estimator_weights(n, "S2", "sd", "blue")
  # the same quartiles as stats::quantile(), type 7, on the first sample
  expect_lt(max(abs(Y[1, ] - quantile(x[1, ], c(0.25, 0.5, 0.75), type = 7, names = FALSE))), 1e-12)
  mu_hat <- drop(Y %*% w_mu); sd_hat <- drop(Y %*% w_sd)
  cm <- c_mu(n, "S2"); cs <- c_sigma(n, "S2")
  # Monte Carlo standard errors: of a mean, SD/sqrt(R); of an SD, about SD/sqrt(2R)
  expect_lt(abs(mean(mu_hat)), 3 * cm / sqrt(R))
  expect_lt(abs(mean(sd_hat) - 1), 3 * cs / sqrt(R))
  expect_lt(abs(sd(mu_hat) - cm), 3 * cm / sqrt(2 * R))
  expect_lt(abs(sd(sd_hat) - cs), 3 * cs / sqrt(2 * R))
})

test_that("the quantile definition is checked, and a rounding residue of the weights adds no rank", {
  expect_true(bad(c_mu(40L, "S2", type = 0)))
  expect_true(bad(c_mu(40L, "S2", type = 10)))
  expect_true(bad(se_imputed(c(q1 = 10, med = 12, q3 = 15), 40L, "S2", type = "7")))
  # under type 8 at n = 41, os_quantile_weights() gives the median a weight of order 1e-15
  # on rank 20 beside 1 on rank 21; the design drops it
  des <- fivenumMeta:::.reported_design(41L, "S2", 8)
  expect_false(20L %in% des$ranks)
  expect_lt(max(abs(colSums(des$W) - 1)), 1e-12)
  # every definition gives a factor, and they differ at an n not of the form 4Q + 1
  v <- vapply(1:9, function(t) c_sigma(40L, "S2", type = t), numeric(1))
  expect_true(all(is.finite(v))); expect_gt(max(v) / min(v) - 1, 1e-3)
})
