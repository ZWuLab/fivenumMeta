# ==============================================================================
# File: test-pool_imputed.R
# Purpose: pool_imputed(): agreement with metafor::rma() under the common-effect
#          and the DerSimonian-Laird model on the two cases of the vignette and a
#          toy set (estimate, standard error, interval, p-value, Q, I^2, tau^2,
#          weights), the closed forms at K = 2, the dropped rows and the errors,
#          the print. The paper's values are in test-vitd.R.

toy <- data.frame(yi = c(0.10, 0.30, -0.20, 0.05),
                  sei_naive = c(0.10, 0.14, 0.12, 0.20),
                  sei = c(0.12, 0.15, 0.13, 0.25))
cols <- c("estimate", "se", "ci.lb", "ci.ub", "pval", "Q", "I2", "tau2")
row_of <- function(f, r) unlist(f$fits[r, cols], use.names = FALSE)

test_that("pool_imputed() agrees with metafor::rma() under EE and DL, both rows, two levels", {
  skip_if_not_installed("metafor")
  cases <- list(toy, se_imputed_diff(vitd, suffixes = c(".controls", ".cases")))
  if (requireNamespace("metamedian", quietly = TRUE))
    cases <- c(cases, list(se_imputed_diff(metamedian::dat.age, measure = "lnVR")))
  for (d in cases) for (m in c("EE", "DL")) for (lev in c(0.95, 0.90)) {
    f <- pool_imputed(d, method = m, level = lev)
    expect_s3_class(f, "pool_imputed")
    expect_identical(f$k, nrow(d)); expect_identical(f$method, m); expect_identical(f$level, lev)
    for (r in c("naive", "exact")) {
      sei <- if (r == "naive") d$sei_naive else d$sei
      ref <- metafor::rma(yi = d$yi, sei = sei, method = m, level = 100 * lev)
      expect_equal(row_of(f, r),
                   c(as.numeric(ref$b), ref$se, ref$ci.lb, ref$ci.ub, ref$pval, ref$QE, ref$I2, ref$tau2),
                   tolerance = 1e-8)
      expect_equal(f$weights[[r]], unname(stats::weights(ref)), tolerance = 1e-8)   # metafor's method, percent
    }
  }
})

test_that("the closed forms hold at K = 2, and the structure of the result", {
  two <- data.frame(yi = c(1, 3), sei_naive = c(0.5, 0.5), sei = c(1, 1))
  f <- pool_imputed(two)
  # equal variances v: the estimate is the mean, the SE sqrt(v / 2), Q = (y1 - y2)^2 / (2 v)
  expect_equal(row_of(f, "naive")[1:2], c(2, sqrt(0.25 / 2)))
  expect_equal(row_of(f, "exact")[1:2], c(2, sqrt(1 / 2)))
  expect_equal(f$fits["naive", "Q"], 4 / (2 * 0.25)); expect_equal(f$fits["exact", "Q"], 4 / 2)
  expect_equal(f$fits["exact", "I2"], 100 * (2 - 1) / 2)
  z <- qnorm(0.975)
  expect_equal(f$fits["exact", "ci.lb"], 2 - z * sqrt(1 / 2)); expect_equal(f$fits["exact", "ci.ub"], 2 + z * sqrt(1 / 2))
  expect_equal(f$fits["exact", "pval"], 2 * pnorm(-2 / sqrt(1 / 2)))
  expect_identical(f$fits$tau2, c(0, 0))                   # the common-effect model sets tau^2 to zero
  expect_identical(rownames(f$fits), c("naive", "exact")); expect_identical(names(f$fits), cols)
  expect_identical(names(f$weights), c("naive", "exact"))
  expect_equal(colSums(f$weights), c(naive = 100, exact = 100))
  expect_true(is.na(f$measure))                             # no attribute on a hand-made data frame
  # DerSimonian-Laird: Q = 2 > K - 1 = 1 under the exact SEs, so tau^2 = (Q - 1) / (sum w - sum w^2 / sum w) = 1
  g <- pool_imputed(two, method = "DL")
  expect_equal(g$fits["exact", "tau2"], (2 - 1) / (2 - 2 / 2))
  expect_equal(g$fits["exact", "se"], sqrt(1 / (2 * (1 / (1 + 1)))))   # two weights 1 / (v + tau^2)
  expect_equal(g$fits["exact", "I2"], f$fits["exact", "I2"])   # Q and I^2 come from the fixed weights
  expect_equal(g$fits["exact", "Q"], f$fits["exact", "Q"])
  # all estimates equal: Q = 0, I^2 = 0, tau^2 = 0
  same <- data.frame(yi = c(1, 1, 1), sei_naive = c(0.2, 0.3, 0.4), sei = c(0.3, 0.4, 0.5))
  h <- pool_imputed(same, method = "DL")
  expect_identical(unname(h$fits$Q), c(0, 0)); expect_identical(unname(h$fits$I2), c(0, 0))
  expect_identical(unname(h$fits$tau2), c(0, 0))
})

test_that("rows that cannot be used are dropped with one warning, and the errors", {
  bad <- rbind(toy, data.frame(yi = c(NA, 0.2, 0.1), sei_naive = c(0.1, NA, 0.1), sei = c(0.1, 0.1, 0)))
  expect_warning(f <- pool_imputed(bad), "3 of 7")
  expect_identical(f$k, 4L)
  expect_equal(f$fits, pool_imputed(toy)$fits)
  expect_identical(nrow(f$weights), 4L)
  expect_error(pool_imputed(toy[1, ]), "at least two")
  expect_error(pool_imputed(toy[, c("yi", "sei")]), "sei_naive")
  expect_error(pool_imputed(as.matrix(toy)), "data frame")
  expect_error(pool_imputed(toy, level = 95)); expect_error(pool_imputed(toy, level = c(0.9, 0.95)))
  expect_error(pool_imputed(toy, method = "REML"))
})

test_that("the print names the measure, the model and the studies, and shows tau2 under DL only", {
  f <- pool_imputed(toy)
  expect_output(print(f), "Pooled estimate")
  expect_output(print(f), "common-effect model")
  expect_output(print(f), "4 studies")
  expect_false(any(grepl("tau2", capture.output(print(f)))))
  g <- pool_imputed(toy, method = "DL")
  expect_output(print(g), "DerSimonian-Laird")
  expect_true(any(grepl("tau2", capture.output(print(g)))))
  d <- se_imputed_diff(vitd, suffixes = c(".controls", ".cases"))
  expect_identical(pool_imputed(d)$measure, "MD")
  expect_output(print(pool_imputed(d)), "difference of means")
  expect_output(print(pool_imputed(d, level = 0.9)), "90%")
  v <- toy; attr(v, "measure") <- "lnVR"
  expect_output(print(pool_imputed(v)), "ratio of SDs")
  out <- capture.output(r <- withVisible(print(f))); expect_false(r$visible)
})

test_that("a factor or character column is refused, and the warning names the dropped rows and their columns", {
  f <- toy; f$yi <- factor(f$yi)
  expect_error(pool_imputed(f), "yi must be numeric, not factor")
  g <- toy; g$sei <- as.character(g$sei)
  expect_error(pool_imputed(g), "sei must be numeric, not character")
  bad <- rbind(toy, data.frame(yi = c(NA, 0.2, 0.1, -0.3), sei_naive = c(0.1, NA, 0.1, 0.1),
                               sei = c(0.1, 0.1, 0, 0.1)))
  expect_warning(h <- pool_imputed(bad), "3 of 8.*5 \\(yi\\), 6 \\(sei_naive\\), 7 \\(sei\\)")
  expect_identical(h$k, 5L)   # a negative effect is a valid effect
})

test_that("the fits are stable when one study dominates, equivariant to the unit of yi, finite next to level = 1", {
  # two studies: tau^2 = max(0, {(y1 - y2)^2 - v1 - v2} / 2), then the weights 1 / (v + tau^2)
  dl2 <- function(y, v) {
    tau2 <- max(0, ((y[1] - y[2])^2 - sum(v)) / 2); w <- 1 / (v + tau2)
    c(sum(w * y) / sum(w), 1 / sqrt(sum(w)), tau2)
  }
  d <- data.frame(yi = c(0, 10), sei_naive = c(1e-8, 1), sei = c(1e-8, 1))
  f <- pool_imputed(d, method = "DL")$fits["exact", ]
  expect_equal(c(f$estimate, f$se, f$tau2), dl2(c(0, 10), c(1e-16, 1)), tolerance = 1e-10)
  base <- data.frame(yi = c(0, 10, 4), sei_naive = c(1, 2, 1.5), sei = c(1.2, 2, 1.5))
  for (m in c("EE", "DL")) {
    f0 <- pool_imputed(base, method = m)$fits
    for (s in c(1e-100, 1e100)) {
      fs <- pool_imputed(base * s, method = m)$fits
      expect_equal(fs[, c("estimate", "se", "ci.lb", "ci.ub")],
                   f0[, c("estimate", "se", "ci.lb", "ci.ub")] * s, tolerance = 1e-10)
      expect_equal(fs$tau2, f0$tau2 * s^2, tolerance = 1e-10)
      expect_equal(fs[, c("pval", "Q", "I2")], f0[, c("pval", "Q", "I2")], tolerance = 1e-10)
    }
  }
  near1 <- pool_imputed(base, level = 1 - .Machine$double.eps / 2)$fits
  expect_true(all(is.finite(unlist(near1))))
  expect_error(pool_imputed(data.frame(yi = c(0, 1e200), sei_naive = c(1, 1), sei = c(1, 1))),
               "not representable")
})

