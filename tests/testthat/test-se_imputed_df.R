# ==============================================================================
# File: test-se_imputed_df.R
# Purpose: se_imputed_df() and se_imputed_diff(): row-by-row agreement with
#          se_imputed(), the sample rows, the mean-with-range rows, the preference
#          switch, the rows that cannot be classified, the suffixes, the difference
#          and the log variability ratio, and the public data set shipped with
#          metamedian.

toy <- data.frame(id = 1:5,
                  n   = c(25, 30, 25, 41, 25),
                  min = c(3, NA, 2, 1, NA), q1 = c(10, 10, NA, 9, NA),
                  med = c(12, 12, 11, 12, NA), q3 = c(15, 15, NA, 16, NA),
                  max = c(30, NA, 25, 30, NA),
                  mean = c(NA, NA, NA, 12.5, 13), sd = c(NA, NA, NA, 4, 3))
# rows: 1 S3 at an exact-rank n; 2 S2 at n = 30 (interpolated quartiles); 3 S1; 4 a summary and
# a sample mean and SD; 5 a sample mean and SD only

test_that("se_imputed_df() reproduces se_imputed() row by row, both versions, and reads the scenarios", {
  for (v in c("blue", "wan_luo_shi")) {
    d <- se_imputed_df(toy, version = v)
    expect_identical(d$scenario, c("S3", "S2", "S1", "sample", "sample"))
    expect_identical(nrow(d), 5L); expect_true(all(c("id", "mu_hat", "exact") %in% names(d)))
    rows <- list(list(1, c(min = 3, q1 = 10, med = 12, q3 = 15, max = 30), 25L, "S3"),
                 list(2, c(q1 = 10, med = 12, q3 = 15), 30L, "S2"),
                 list(3, c(min = 2, med = 11, max = 25), 25L, "S1"))
    for (r in rows) {
      i <- r[[1]]; ref <- se_imputed(r[[2]], r[[3]], r[[4]], version = v)
      for (f in c("mu_hat", "sigma_hat", "se_naive", "se_mu", "se_sigma"))
        expect_lt(abs(d[[f]][i] - ref[[f]]), 1e-12)
      expect_identical(d$version_used[i], ref$version_used)
      expect_identical(d$exact[i], ref$exact)
    }
    expect_true(isTRUE(d$exact[2])); expect_identical(d$version_used[2], v)
  }
})

test_that("a sample row gives the sample mean's SE and the SD of the sample SD, and prefer switches it", {
  d <- se_imputed_df(toy)
  bS <- function(n) sqrt(2 / (n - 1)) * exp(lgamma(n / 2) - lgamma((n - 1) / 2))
  for (i in 4:5) {
    expect_identical(d$version_used[i], "sample"); expect_true(is.na(d$exact[i]))
    expect_identical(d$mu_hat[i], toy$mean[i]); expect_identical(d$sigma_hat[i], toy$sd[i])
    expect_lt(abs(d$se_mu[i] - toy$sd[i] / sqrt(toy$n[i])), 1e-14)
    expect_lt(abs(d$se_naive[i] - toy$sd[i] / sqrt(toy$n[i])), 1e-14)
    expect_lt(abs(d$se_sigma[i] - toy$sd[i] * sqrt(1 - bS(toy$n[i])^2)), 1e-12)
  }
  # prefer = "imputed": row 4 is imputed from its five numbers, row 5 stays a sample row
  d2 <- se_imputed_df(toy, prefer = "imputed")
  expect_identical(d2$scenario, c("S3", "S2", "S1", "S3", "sample"))
  expect_error(se_imputed_df(toy, prefer = "reported"))   # the former name of the default is gone
  ref <- se_imputed(c(min = 1, q1 = 9, med = 12, q3 = 16, max = 30), 41L, "S3")
  expect_lt(abs(d2$mu_hat[4] - ref$mu_hat), 1e-12); expect_true(isTRUE(d2$exact[4]))
})

test_that("rows that cannot be used get NA and one warning, and the function does not stop", {
  bad <- data.frame(n = c(25, NA, 25, 25.5),
                    q1 = c(10, 10, 15, 10), med = c(12, 12, 12, 12), q3 = c(15, 15, 10, 15))
  # row 2: no sample size; row 3: values out of order; row 4: a non-integer n
  expect_warning(d <- se_imputed_df(bad), "3 of 4")
  expect_false(is.na(d$mu_hat[1])); expect_true(all(is.na(d$mu_hat[2:4])))
  expect_identical(d$scenario[1], "S2"); expect_true(all(is.na(d$scenario[2:4])))
  # a row with nothing usable
  none <- data.frame(n = 25, med = 12)
  expect_warning(d3 <- se_imputed_df(none), "no usable summary")
  expect_true(is.na(d3$scenario))
  # not a data frame, a non-numeric column
  expect_error(se_imputed_df(c(n = 25)))
  expect_error(se_imputed_df(data.frame(n = 25, q1 = "10", med = 12, q3 = 15)))
})

test_that("suffixes select the columns and name the outputs, absent columns count as absent, append = FALSE returns the outputs alone", {
  two <- data.frame(n.g1 = c(25, 30), q1.g1 = c(10, 11), med.g1 = c(12, 13), q3.g1 = c(15, 16),
                    n.g2 = c(29, 40), q1.g2 = c(8, 9), med.g2 = c(10, 11), q3.g2 = c(13, 14))
  d1 <- se_imputed_df(two, suffix = ".g1", append = FALSE)
  expect_identical(names(d1), paste0(c("scenario", "mu_hat", "sigma_hat", "se_naive", "se_mu",
                                       "se_sigma", "version_used", "exact"), ".g1"))
  expect_identical(d1$scenario.g1, c("S2", "S2"))      # no min/max/mean/sd columns: S2 from the quartiles
  ref <- se_imputed(c(q1 = 10, med = 12, q3 = 15), 25L, "S2")
  expect_lt(abs(d1$mu_hat.g1[1] - ref$mu_hat), 1e-12)
  # a column argument set to NULL is ignored even when the column exists
  d4 <- se_imputed_df(toy[1, ], min = NULL, max = NULL)
  expect_identical(d4$scenario, "S2")
})

test_that("se_imputed_diff() gives group 1 minus group 2 with the variances added", {
  two <- data.frame(n.g1 = c(25, 30), q1.g1 = c(10, 11), med.g1 = c(12, 13), q3.g1 = c(15, 16),
                    n.g2 = c(29, 40), q1.g2 = c(8, 9), med.g2 = c(10, 11), q3.g2 = c(13, 14))
  d <- se_imputed_diff(two)
  g1 <- se_imputed_df(two, suffix = ".g1", append = FALSE)
  g2 <- se_imputed_df(two, suffix = ".g2", append = FALSE)
  expect_lt(max(abs(d$yi - (g1$mu_hat.g1 - g2$mu_hat.g2))), 1e-12)
  expect_lt(max(abs(d$sei - sqrt(g1$se_mu.g1^2 + g2$se_mu.g2^2))), 1e-12)
  expect_lt(max(abs(d$sei_naive - sqrt(g1$se_naive.g1^2 + g2$se_naive.g2^2))), 1e-12)
  expect_true(all(d$sei > d$sei_naive))         # the exact SE exceeds the naive one under S2
  expect_error(se_imputed_diff(two, suffixes = ".g1"))
  expect_error(se_imputed_diff(two, suffixes = c(".g1", ".g1")))
})

test_that("a mean with the range and no SD is the scenario mean_range: the mean kept, the SD from the range", {
  mr <- data.frame(n = c(35, 16, 25, 25, 25), mean = c(26.75, 48.5, 13, 13, 13),
                   sd = c(NA, NA, 3, NA, NA), min = c(2.5, 22.5, 2, 2, 5),
                   max = c(75, 145, 25, 25, 5), med = c(NA, NA, NA, 11, NA))
  # rows: 1, 2 the vitamin D study with a mean and the range (cases, controls); 3 a mean, the range
  # and an SD, so a sample row; 4 a mean beside a complete S1 summary, so S1 with the mean ignored;
  # 5 a zero range
  for (v in c("blue", "wan_luo_shi")) {
    expect_warning(d <- se_imputed_df(mr, version = v), "zero range")
    expect_identical(d$scenario, c("mean_range", "mean_range", "sample", "S1", NA))
    for (i in 1:2) {
      # the S1 SD estimator with the midrange as the median: its weight on the median is zero
      ref <- se_imputed(c(min = mr$min[i], med = (mr$min[i] + mr$max[i]) / 2, max = mr$max[i]),
                        mr$n[i], "S1", version = v)
      expect_lt(abs(ref$w_sd[2]), 1e-12)
      expect_identical(d$mu_hat[i], mr$mean[i])
      expect_lt(abs(d$sigma_hat[i] - ref$sigma_hat), 1e-12)
      expect_lt(abs(d$se_naive[i] - ref$sigma_hat / sqrt(mr$n[i])), 1e-14)
      expect_lt(abs(d$se_mu[i] - ref$sigma_hat / sqrt(mr$n[i])), 1e-14)
      expect_lt(abs(d$se_sigma[i] - ref$se_sigma), 1e-12)
      expect_identical(d$version_used[i], v); expect_true(isTRUE(d$exact[i]))
    }
    # the BLUE: the range over the expected range of a normal sample, the paper's xi(n)
    if (v == "blue") expect_equal(d$sigma_hat[1:2], c(17.207745925042, 34.6830682419412), tolerance = 1e-6)
    expect_identical(d$version_used[3], "sample"); expect_identical(d$sigma_hat[3], 3)
    ref4 <- se_imputed(c(min = 2, med = 11, max = 25), 25L, "S1", version = v)
    expect_lt(abs(d$mu_hat[4] - ref4$mu_hat), 1e-12)
    expect_true(is.na(d$scenario[5]) && is.na(d$mu_hat[5]) && is.na(d$se_sigma[5]))
  }
  # prefer = "imputed" changes nothing for a mean_range row: nothing is imputed to prefer
  expect_identical(se_imputed_df(mr[1:2, ], prefer = "imputed")$scenario, c("mean_range", "mean_range"))
  # a non-finite mean, and a mean with the minimum alone (no usable summary)
  expect_warning(d6 <- se_imputed_df(data.frame(n = c(25, 25), mean = c(Inf, 10), min = c(1, 1), max = c(9, NA))), "2 of 2")
  expect_true(all(is.na(d6$scenario)))
})

test_that("measure = \"lnVR\" gives the log ratio of the SDs with the two variances of the paper, and the attribute", {
  two <- data.frame(n.g1 = c(25, 30, 20), q1.g1 = c(10, 11, NA), med.g1 = c(12, 13, NA),
                    q3.g1 = c(15, 16, NA), mean.g1 = c(NA, NA, 10), sd.g1 = c(NA, NA, 3),
                    n.g2 = c(29, 40, 22), q1.g2 = c(8, 9, 7), med.g2 = c(10, 11, 9), q3.g2 = c(13, 14, 12))
  d <- se_imputed_diff(two, measure = "lnVR"); m <- se_imputed_diff(two)
  expect_identical(attr(d, "measure"), "lnVR"); expect_identical(attr(m, "measure"), "MD")
  expect_equal(d$yi, log(d$sigma_hat.g1 / d$sigma_hat.g2))
  expect_equal(d$sei_naive^2, 1 / (2 * (two$n.g1 - 1)) + 1 / (2 * (two$n.g2 - 1)))
  # group 1: two imputed terms C_sigma(n)^2 and one sample term 1 / (2 (n - 1)); group 2: imputed
  t1 <- ifelse(d$scenario.g1 == "sample", 1 / (2 * (two$n.g1 - 1)), (d$se_sigma.g1 / d$sigma_hat.g1)^2)
  t2 <- (d$se_sigma.g2 / d$sigma_hat.g2)^2
  expect_equal(d$sei^2, t1 + t2)
  expect_equal(t2[1], c_sigma(29, "S2")^2, tolerance = 1e-10)
  expect_true(all(d$sei > d$sei_naive))
  # the per-group columns are those of the default measure
  for (cl in setdiff(names(d), c("yi", "sei_naive", "sei"))) expect_identical(d[[cl]], m[[cl]])
  # the sample sizes are read from the column the argument n names
  names(two) <- sub("^n[.]", "N.", names(two))
  d2 <- se_imputed_diff(two, n = "N", measure = "lnVR")
  expect_equal(d2$sei, d$sei); expect_equal(d2$sei_naive, d$sei_naive)
  suppressWarnings(expect_error(se_imputed_diff(two, measure = "lnVR"), "sample size"))
  expect_error(se_imputed_diff(two, n = "N", measure = "VR"))
})

test_that("the public data set of metamedian goes through: every row classified, no NA, and rma() runs", {
  skip_if_not_installed("metamedian")
  dat <- metamedian::dat.age
  d <- se_imputed_diff(dat)
  expect_identical(nrow(d), nrow(dat))
  expect_identical(as.integer(table(d$scenario.g1)[c("S2", "sample")]), c(29L, 22L))
  expect_false(anyNA(d$yi)); expect_false(anyNA(d$sei)); expect_false(anyNA(d$sei_naive))
  expect_true(all(d$sei >= d$sei_naive))
  skip_if_not_installed("metafor")
  fit <- metafor::rma(yi, sei = sei, data = d)
  expect_s3_class(fit, "rma")
})

test_that("a rerun on a data frame that carries the output columns replaces them, and se_imputed_diff() reads the new ones", {
  two <- data.frame(n.g1 = 25, q1.g1 = 10, med.g1 = 12, q3.g1 = 15, n.g2 = 29, q1.g2 = 8, med.g2 = 10, q3.g2 = 13)
  old <- se_imputed_diff(two)
  old$med.g1 <- 14; two$med.g1 <- 14
  rerun <- se_imputed_diff(old); fresh <- se_imputed_diff(two)
  expect_false(anyDuplicated(names(rerun)) > 0)
  expect_identical(rerun$yi, fresh$yi); expect_identical(rerun$sei, fresh$sei)
  # a change of version, type or prefer on the earlier result equals the fresh call
  for (args in list(list(version = "wan_luo_shi"), list(type = 6), list(prefer = "imputed"))) {
    a <- do.call(se_imputed_diff, c(list(old), args)); b <- do.call(se_imputed_diff, c(list(two), args))
    expect_identical(a[names(b)], b[names(b)])   # the same columns; [ drops the measure attribute of both
  }
  d1 <- se_imputed_df(toy); d2 <- se_imputed_df(d1)
  expect_identical(names(d2), names(d1)); expect_identical(d2, d1)
})

test_that("non-finite reported rows and sample sizes get NA and the one warning, with no coercion warning", {
  n_warn <- function(expr) { k <- 0L; withCallingHandlers(expr, warning = function(w) { k <<- k + 1L; invokeRestart("muffleWarning") }); k }
  d <- data.frame(n = c(25, 25, 25, Inf, 2^31), mean = c(10, Inf, 10, 10, 10), sd = c(3, 3, Inf, 3, 3))
  expect_identical(n_warn(r <- se_imputed_df(d, append = FALSE)), 1L)
  expect_warning(se_imputed_df(d), "4 of 5")
  expect_false(anyNA(r[1, c("mu_hat", "sigma_hat", "se_naive", "se_mu", "se_sigma")]))   # exact is NA for a reported row
  expect_true(all(is.na(r$mu_hat[2:5]))); expect_true(all(is.na(r$scenario[2:5])))
  expect_warning(se_imputed_df(data.frame(n = Inf, q1 = 10, med = 12, q3 = 15)), "sample size")
})

test_that("a zero-row data frame comes back with the output schema, and an all-NA logical column counts as absent", {
  empty <- data.frame(n = numeric(), q1 = numeric(), med = numeric(), q3 = numeric())
  d0 <- se_imputed_df(empty)
  expect_identical(nrow(d0), 0L)
  expect_identical(names(d0), c(names(empty), "scenario", "mu_hat", "sigma_hat", "se_naive", "se_mu", "se_sigma", "version_used", "exact"))
  expect_identical(nrow(se_imputed_df(empty, append = FALSE)), 0L)
  two0 <- data.frame(n.g1 = numeric(), q1.g1 = numeric(), med.g1 = numeric(), q3.g1 = numeric(),
                     n.g2 = numeric(), q1.g2 = numeric(), med.g2 = numeric(), q3.g2 = numeric())
  expect_identical(nrow(se_imputed_diff(two0)), 0L)
  # read.csv() types an all-NA column as logical
  d <- data.frame(n = c(25, 30), q1 = c(10, 11), med = c(12, 13), q3 = c(15, 16), min = NA, max = NA)
  expect_identical(class(d$min), "logical")
  r <- se_imputed_df(d, append = FALSE)
  expect_identical(r$scenario, c("S2", "S2"))
  expect_lt(abs(r$mu_hat[1] - se_imputed(c(q1 = 10, med = 12, q3 = 15), 25L, "S2")$mu_hat), 1e-12)
})

test_that("a reported SD of zero and a constant summary get NA with their reason, in the one warning", {
  n_warn <- function(expr) { k <- 0L; withCallingHandlers(expr, warning = function(w) { k <<- k + 1L; invokeRestart("muffleWarning") }); k }
  d <- data.frame(n = c(25, 25, 25, 30),
                  q1 = c(10, 5, NA, 7), med = c(12, 5, NA, 7), q3 = c(15, 5, NA, 7),
                  min = c(NA, NA, NA, 7), max = c(NA, NA, NA, 7),
                  mean = c(NA, NA, 10, NA), sd = c(NA, NA, 0, NA))
  # rows: 1 S2; 2 S2 with all three quartiles equal; 3 a reported pair with sd = 0; 4 S3, all five equal
  expect_identical(n_warn(r <- se_imputed_df(d, append = FALSE)), 1L)   # no per-row warning escapes
  expect_warning(se_imputed_df(d), "3 of 4: 2 (constant summary); 3 (zero SD); 4 (constant summary)", fixed = TRUE)
  expect_false(anyNA(r[1, c("mu_hat", "sigma_hat", "se_mu", "se_sigma")]))
  expect_true(all(is.na(r$mu_hat[2:4]))); expect_true(all(is.na(r$se_mu[2:4]))); expect_true(all(is.na(r$scenario[2:4])))
  # the two-group function carries the NA into the difference
  two <- data.frame(n.g1 = 25, q1.g1 = 5, med.g1 = 5, q3.g1 = 5, n.g2 = 25, q1.g2 = 10, med.g2 = 12, q3.g2 = 15)
  expect_warning(dd <- se_imputed_diff(two), "constant summary")
  expect_true(is.na(dd$yi)); expect_true(is.na(dd$sei)); expect_false(is.na(dd$mu_hat.g2))
})

test_that("a mean_range row uses the two extremes alone: one SD under every quantile type, even and odd n, against the extremes' moments", {
  for (n in c(10L, 11L, 16L, 35L)) for (v in c("blue", "wan_luo_shi")) {
    mo <- osmoments::os_moments(n, c(1L, n))
    div <- unname(if (v == "blue") mo$alpha[2] - mo$alpha[1] else 2 * qnorm((n - 0.375) / (n + 0.25)))
    sd_range <- sqrt(mo$B[1, 1] + mo$B[2, 2] - 2 * mo$B[1, 2])
    x <- data.frame(n = n, min = 1, max = 9, mean = 5)
    for (t in 1:9) {
      r <- se_imputed_df(x, version = v, type = t, append = FALSE)
      expect_equal(r$sigma_hat, 8 / div, tolerance = 1e-12)
      expect_equal(r$se_sigma, 8 / div * sd_range / div, tolerance = 1e-12)
      expect_equal(r$se_mu, r$sigma_hat / sqrt(n)); expect_equal(r$se_naive, r$se_mu)
      expect_identical(r$version_used, v); expect_true(r$exact)
    }
  }
  # the rows refused, each with its reason, in one warning
  bad <- data.frame(n = c(25, 25, 25, 2, 25), min = c(9, 1, 1, 1, 1), max = c(1, 1, 9, 9, 9),
                    mean = c(5, 1, 50, 5, -2))
  expect_warning(d <- se_imputed_df(bad),
                 "1 \\(minimum above maximum\\); 2 \\(zero range\\); 3 \\(mean outside the range\\); 4 \\(a mean with the range needs n >= 3\\); 5 \\(mean outside the range\\)")
  expect_true(all(is.na(d$scenario)))
})

test_that("under lnVR the sample sizes follow the column mapping of se_imputed_df(), and an unusable group gives NA in all three columns", {
  d <- data.frame(N.g1 = 25, mean.g1 = 1, sd.g1 = 2, N.g2 = 30, mean.g2 = 1, sd.g2 = 1,
                  n.g1 = 1000, n.g2 = 1000)   # competing columns under the default name
  ref <- sqrt(1 / (2 * 24) + 1 / (2 * 29))
  by_name <- se_imputed_diff(d, n = "N", measure = "lnVR")
  by_pos  <- se_imputed_diff(d, c(".g1", ".g2"), "N", measure = "lnVR")
  expect_equal(by_name$yi, log(2)); expect_equal(by_name$sei_naive, ref); expect_equal(by_name$sei, ref)
  expect_identical(by_pos[, c("yi", "sei_naive", "sei")], by_name[, c("yi", "sei_naive", "sei")])
  suppressWarnings(expect_error(se_imputed_diff(d, n = NULL, measure = "lnVR"), "sample size"))
  # a zero SD or an invalid sample size in one group: NA in yi, sei_naive and sei alike
  z <- data.frame(n.g1 = c(2, 0, 1, 1.5, 10), mean.g1 = 1, sd.g1 = c(0, 2, 2, 2, 1e300),
                  n.g2 = c(3, 25, 25, 25, 10), mean.g2 = 1, sd.g2 = c(2, 1, 1, 1, 1e-300))
  expect_warning(r <- se_imputed_diff(z, measure = "lnVR"), "4 of 5")
  expect_true(all(is.na(r$yi[1:4])), all(is.na(r$sei_naive[1:4])), all(is.na(r$sei[1:4])))
  # the log ratio as a difference of logs, finite where the ratio itself overflows
  expect_equal(r$yi[5], log(1e300) - log(1e-300)); expect_true(is.finite(r$sei[5]))
})

