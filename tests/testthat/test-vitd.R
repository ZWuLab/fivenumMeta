# ==============================================================================
# File: test-vitd.R
# Purpose: The data set vitd as typed from its source, and the two meta-analyses of
#          the illustration section of Wu, Shi and Tong (2026) (its Tables 5 and 6)
#          reproduced with the package in three lines each: the mean case on vitd
#          and the SD case on metamedian::dat.age, each pooled under the naive and
#          the exact standard errors with the common-effect model by pool_imputed().
#          The expected values are the full-precision values that the analysis
#          behind that section wrote to its output files; the paper's printed values
#          follow by rounding.

tol <- 1e-6
fit_row <- function(f, r, cols) unlist(f$fits[r, cols], use.names = FALSE)

test_that("vitd is the six studies as typed from Table 1 of Nnoaham and Clarke (2008)", {
  expect_s3_class(vitd, "data.frame")
  expect_equal(dim(vitd), c(6L, 13L))
  expect_equal(names(vitd),
               c("study", paste0(c("n", "min", "med", "max", "mean", "sd"), ".cases"),
                 paste0(c("n", "min", "med", "max", "mean", "sd"), ".controls")))
  expect_equal(vitd$n.cases, c(40, 40, 15, 51, 22, 35))
  expect_equal(vitd$n.controls, c(40, 38, 15, 51, 23, 16))
  # the three studies that report the median and the range, and no mean or SD
  s1 <- !is.na(vitd$med.cases)
  expect_equal(which(s1), 1:3)
  expect_true(all(is.na(vitd$mean.cases[s1])) && all(is.na(vitd$sd.cases[s1])))
  # the sums of the reported values, a check of the typing
  expect_equal(sum(vitd[1:3, c("min.cases", "med.cases", "max.cases")]), 478.25)
  expect_equal(sum(vitd[1:3, c("min.controls", "med.controls", "max.controls")]), 618.25)
  expect_equal(sum(vitd$mean.cases, na.rm = TRUE), 142.75)
  expect_equal(sum(vitd$mean.controls, na.rm = TRUE), 196.25)
  expect_equal(sum(vitd$sd.cases, na.rm = TRUE), 43)
  expect_equal(sum(vitd$sd.controls, na.rm = TRUE), 45)
  # the study with a mean and the range has no SD and no median
  i <- vitd$study == "Sasidharan, India"
  expect_true(is.na(vitd$sd.cases[i]) && is.na(vitd$med.cases[i]) && !is.na(vitd$min.cases[i]))
})

test_that("the mean case reproduces Table 5 of the paper, with the data set as the source prints it", {
  # controls minus cases, the paper's difference
  d <- se_imputed_diff(vitd, suffixes = c(".controls", ".cases"))
  expect_identical(attr(d, "measure"), "MD")
  expect_equal(d$scenario.cases, c("S1", "S1", "S1", "sample", "sample", "mean_range"))
  expect_equal(d$scenario.controls, d$scenario.cases)
  expect_true(all(d$exact.cases[1:3]))                       # the exact BLUE at n = 40, 40, 15
  # the study with a mean and the range: the SD is the range over xi(n), the mean is kept
  i <- which(d$study == "Sasidharan, India")
  expect_equal(d$sigma_hat.cases[i], 17.207745925042, tolerance = tol)
  expect_equal(d$sigma_hat.controls[i], 34.6830682419412, tolerance = tol)
  expect_identical(d$mu_hat.cases[i], vitd$mean.cases[i]); expect_identical(d$mu_hat.controls[i], vitd$mean.controls[i])
  expect_true(isTRUE(d$exact.cases[i])); expect_identical(d$version_used.cases[i], "blue")
  expect_equal(d$yi, c(15.4246439635143, 3.03986754287223, 22.9126394532552, 26, 5.75, 21.75),
               tolerance = tol)
  expect_equal(d$sei_naive, c(5.23034850982935, 4.29888594702797, 8.50826967370674,
                              5.34278421211777, 5.1324598652863, 9.14562123082905), tolerance = tol)
  expect_equal(d$sei, c(6.00637836197536, 4.93055145295272, 9.33953705489607,
                        5.34278421211777, 5.1324598652863, 9.14562123082905), tolerance = tol)
  # the common-effect fit under each standard error: the paper prints 12.8 (2.30) and 13.3 (2.46)
  f <- pool_imputed(d)
  expect_identical(f$measure, "MD"); expect_identical(f$k, 6L)
  expect_equal(fit_row(f, "naive", c("estimate", "se")), c(12.849079839385, 2.29761112499946), tolerance = tol)
  expect_equal(fit_row(f, "exact", c("estimate", "se")), c(13.338517395062, 2.45983630013422), tolerance = tol)
  expect_equal(round(100 * (f$fits["exact", "se"] / f$fits["naive", "se"] - 1), 1), 7.1)   # the rise the paper quotes
  # the share of the weight on the three S1 studies falls, as the paper reads it
  expect_lt(sum(f$weights$exact[1:3]), sum(f$weights$naive[1:3]))
})

test_that("the SD case reproduces Table 6 of the paper through measure = \"lnVR\"", {
  skip_if_not_installed("metamedian")
  d <- se_imputed_diff(metamedian::dat.age, measure = "lnVR")
  expect_identical(attr(d, "measure"), "lnVR")
  expect_equal(nrow(d), 51L)
  expect_equal(sum(d$scenario.g1 == "S2"), 29L)
  expect_equal(sum(d$scenario.g1 == "sample"), 22L)
  # the effect and its two variances as the paper defines them
  expect_equal(d$yi, log(d$sigma_hat.g1 / d$sigma_hat.g2))
  expect_equal(d$sei_naive^2, 1 / (2 * (d$n.g1 - 1)) + 1 / (2 * (d$n.g2 - 1)))
  term <- function(scenario, n, se_sigma, sigma_hat)
    ifelse(scenario == "sample", 1 / (2 * (n - 1)), (se_sigma / sigma_hat)^2)
  expect_equal(d$sei^2, term(d$scenario.g1, d$n.g1, d$se_sigma.g1, d$sigma_hat.g1) +
                        term(d$scenario.g2, d$n.g2, d$se_sigma.g2, d$sigma_hat.g2))
  # an imputed group's term is C_sigma(n)^2 at its own sample size
  imp <- which(d$scenario.g1 == "S2")[1:3]
  expect_equal((d$se_sigma.g1[imp] / d$sigma_hat.g1[imp])^2,
               sapply(d$n.g1[imp], function(n) c_sigma(n, "S2"))^2, tolerance = 1e-10)
  f <- pool_imputed(d)
  # the paper prints -0.225 (0.0161), Q 346.6, I2 86 and -0.215 (0.0217), Q 165.3, I2 70
  expect_equal(fit_row(f, "naive", c("estimate", "se", "Q", "I2")),
               c(-0.225278160713769, 0.0161016004846129, 346.635900749338, 85.5756429464136), tolerance = tol)
  expect_equal(fit_row(f, "exact", c("estimate", "se", "Q", "I2")),
               c(-0.215248143041705, 0.0216661887161565, 165.275547047304, 69.7474908458845), tolerance = tol)
})
