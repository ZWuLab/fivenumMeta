# ==============================================================================
# File: helper-fivenumMeta.R
# Purpose: Helpers shared by the test files, from the former test script
#          R/se_imputed_test.R of the research repository (commit db76e19).

`%~%` <- function(x, y) max(abs(x - y))   # max absolute difference
bad <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
warned <- function(expr) inherits(tryCatch(expr, warning = function(w) w), "warning")
meth_mu <- c(blue = "blue", wan_luo_shi = "luo")
meth_sd <- function(v, sc) if (v == "blue") "blue" else if (sc == "S3") "shi" else "wan"
