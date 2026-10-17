# ==============================================================================
# File: se_imputed.R
# Purpose: The study-level function: a study's reported values and sample size in, the
#          imputed mean and SD with their plug-in standard errors out.
#
# EXPORTED FUNCTIONS:
#   - se_imputed()        the study-level wrapper, version "blue" (default) or "wan_luo_shi"
#   - print.se_imputed()  the compact print of its result: the estimates, the standard errors, the factors
#
# INTERNAL HELPERS:
#   none

#################### EXPORTED MAIN FUNCTIONS ####################

#' Imputed mean and SD of one study, with their plug-in standard errors
#'
#' The study-level wrapper for practitioners: give it a study's reported values
#' and sample size, and it returns the imputed mean \eqn{\hat\mu} and SD
#' \eqn{\hat\sigma}, the naive standard error \eqn{\hat\sigma/\sqrt n} that
#' treats \eqn{\hat\mu} as a sample mean, and the plug-in standard errors
#' \eqn{\mathrm{SE}(\hat\mu) = C_\mu(n)\,\hat\sigma} and
#' \eqn{\mathrm{SE}(\hat\sigma) = C_\sigma(n)\,\hat\sigma} built from the
#' normal-model factors of \code{\link{c_mu}} and \code{\link{c_sigma}}. The
#' default estimator version is the BLUE (Lloyd's least-squares estimators from
#' the reported values, the exact-optimal version of the Wan-Luo-Shi
#' estimators); \code{version = "wan_luo_shi"} gives the estimators of Wan, Luo
#' and Shi (Luo's mean, Wan's SD under S1 and S2, Shi's under S3, with their
#' closed-form weights and the Blom formulas) and their own factors. Both are
#' computed exactly at every sample size.
#'
#' @param q Named numeric vector of the reported values: \code{min}, \code{med},
#'   \code{max} for S1; \code{q1}, \code{med}, \code{q3} for S2; all five
#'   (\code{min}, \code{q1}, \code{med}, \code{q3}, \code{max}) for S3;
#'   \code{p10}, \code{q1}, \code{med}, \code{q3}, \code{p90} for S4; and
#'   \code{min}, \code{p10}, \code{q1}, \code{med}, \code{q3}, \code{p90},
#'   \code{max} for S5, where \code{p10} and \code{p90} are the 10th and 90th
#'   percentiles. Extra names are ignored, a required name given twice is an
#'   error, and the entries used must be non-decreasing.
#' @param n Integer sample size: at least 3 under S1 and S2, 5 under S3 and S4,
#'   and 7 under S5.
#' @param scenario One of "S1" to "S5".
#' @param version \code{"blue"} (default) or \code{"wan_luo_shi"}.
#' @param type The definition of the reported sample quantiles, an integer from
#'   1 to 9 as in \code{\link[stats]{quantile}}; 7, R's default, unless the
#'   study states another. See Details.
#' @param x An object returned by \code{se_imputed()}.
#' @param digits The number of significant digits the print shows.
#' @param ... Passed to \code{\link{print.default}} for the estimates.
#' @return An object of class \code{"se_imputed"}: a list, used as any list, whose
#'   print shows the estimates, the standard errors and the factors. Its
#'   components: \code{mu_hat}, \code{sigma_hat}, \code{se_naive}
#'   \eqn{= \hat\sigma/\sqrt n}, \code{se_mu} \eqn{= C_\mu(n)\hat\sigma},
#'   \code{se_sigma} \eqn{= C_\sigma(n)\hat\sigma}, the factors \code{c_n} and
#'   \code{c_sigma} used, \code{exact} (TRUE: the factors are computed exactly,
#'   at every \eqn{n}), \code{version} and \code{version_used} (the same, kept
#'   for the columns of \code{\link{se_imputed_df}}), plus \code{n},
#'   \code{scenario}, \code{type} and the weights \code{w_mean}, \code{w_sd},
#'   named by the reported values as \code{\link{estimator_weights}} names them
#'   (\code{a} the minimum, \code{m} the median, \code{b} the maximum,
#'   \code{q1}, \code{q3}, \code{p10}, \code{p90}). The example shows how to
#'   pick the fields a meta-analysis needs.
#' @details
#' The sample size. At \eqn{n = 4Q + 1} with integer \eqn{Q \ge 1}, under the
#' default definition (type 7), the reported quartiles are the order statistics
#' at ranks \eqn{Q+1, 2Q+1, 3Q+1}. At any other \eqn{n}, each reported sample
#' quantile is a fixed linear combination of two adjacent order statistics, and
#' the estimators and their factors are computed exactly from the means and
#' covariance matrix of the reported values (see \code{\link{imputed_moments}}).
#' At such \eqn{n} the factors depend on the definition of the sample quantiles;
#' the median of an even sample is the average of the two middle values under
#' type 7 (and types 2, 5, 6, 8 and 9).
#'
#' The shortcut formulas for the factors (\code{\link{c_mu_shortcut}},
#' \code{\link{c_sigma_shortcut}}) are not used here. They serve a reader
#' without this package and a comparison with the exact values.
#'
#' Under S4 and S5 only the BLUE is defined; \code{version = "wan_luo_shi"} is
#' an error there.
#'
#' A constant summary, all reported values equal, has an imputed SD of exactly
#' zero (the SD weights sum to zero) and so standard errors of zero, which cannot
#' enter inverse-variance weights. The function returns them with a warning;
#' \code{\link{se_imputed_df}} sets such a row to \code{NA} instead.
#'
#' Assumptions: iid normal data. The standard errors are plug-in quantities,
#' \eqn{\sigma} being replaced by \eqn{\hat\sigma}; their reference distribution
#' is not Student's t and is not supplied here.
#' @examples
#' r <- se_imputed(c(q1 = 10, med = 12, q3 = 15), n = 25, scenario = "S2")
#' r                                        # the compact print
#' r[c("mu_hat", "se_mu")]                  # the fields a meta-analysis needs, as a list
#' se_imputed(c(q1 = 10, med = 12, q3 = 15), n = 25, scenario = "S2", version = "wan_luo_shi")
#' se_imputed(c(q1 = 10, med = 12, q3 = 15), n = 24, scenario = "S2")   # interpolated quartiles
#' se_imputed(c(min = 3, q1 = 10, med = 12, q3 = 15, max = 30), n = 30, "S3")
#' @export
se_imputed <- function(q, n, scenario = c("S3", "S1", "S2", "S4", "S5"), version = c("blue", "wan_luo_shi"),
                       type = 7) {
  scenario <- match.arg(scenario); version <- match.arg(version)
  # n: one finite integer >= 2 within R's integer range
  .check_n_integer(n)
  if (n < 2) stop("n must be a single finite integer >= 2; got: ", n)
  n <- as.integer(n)
  need <- switch(scenario, S1 = c("min", "med", "max"),
                 S2 = c("q1", "med", "q3"),
                 S3 = c("min", "q1", "med", "q3", "max"),
                 S4 = c("p10", "q1", "med", "q3", "p90"),
                 S5 = c("min", "p10", "q1", "med", "q3", "p90", "max"))
  # q must be numeric before the names are read: a factor would otherwise be
  # coerced to its level codes by as.numeric() and go through unnoticed
  if (!is.numeric(q) || is.factor(q))
    stop("q must be a named numeric vector, not a ", class(q)[1])
  if (is.null(names(q)) || !all(need %in% names(q)))
    stop("q must be a named numeric vector with entries ",
         paste(need, collapse = ", "), " for scenario ", scenario)
  # a required name given twice is ambiguous; q[need] would take the first silently
  dup <- unique(names(q)[duplicated(names(q)) & names(q) %in% need])
  if (length(dup))
    stop("q names a reported value more than once: ", paste(dup, collapse = ", "))
  x <- as.numeric(q[need])
  if (any(!is.finite(x))) stop("reported values must be finite numbers")
  if (is.unsorted(x))
    stop("reported values must be non-decreasing in the order ",
         paste(need, collapse = " <= "))
  .check_scenario_n(n, scenario)
  if (version == "wan_luo_shi" && scenario %in% c("S4", "S5"))
    stop("the Wan-Luo-Shi estimators are not defined for ", scenario, "; use version = \"blue\"")
  # the exact moments of the reported values, under quantile definition `type`, at any n
  mo <- .reported_moments(n, scenario, type)$moments
  m_mean <- if (version == "blue") "blue" else "luo"
  m_sd   <- if (version == "blue") "blue" else if (scenario == "S3") "shi" else "wan"
  w_mean <- estimator_weights(n, scenario, "mean", m_mean, moments = mo, type = type)
  w_sd   <- estimator_weights(n, scenario, "sd", m_sd, moments = mo, type = type)
  cn <- c_mu(n, scenario, method = m_mean, moments = mo, type = type)
  cs <- c_sigma(n, scenario, method = m_sd, moments = mo, type = type)
  mu_hat <- sum(w_mean * x)
  sigma_hat <- sum(w_sd * x)
  # a constant summary: the SD weights sum to zero, so the imputed SD is zero by
  # construction (not the rounding residue of the sum), and so are the standard errors
  if (x[1] == x[length(x)]) {
    sigma_hat <- 0
    warning("the reported values are all equal, so the imputed SD and the standard errors are zero")
  }
  structure(list(mu_hat = mu_hat, sigma_hat = sigma_hat,
                 se_naive = sigma_hat / sqrt(n),
                 se_mu = cn * sigma_hat, se_sigma = cs * sigma_hat,
                 c_n = cn, c_sigma = cs, exact = TRUE,
                 version = version, version_used = version,
                 n = n, scenario = scenario, type = type, w_mean = w_mean, w_sd = w_sd),
            class = c("se_imputed", "list"))
}

#' @rdname se_imputed
#' @export
print.se_imputed <- function(x, digits = 4L, ...) {
  version <- switch(x$version_used, blue = "the BLUE",
                    wan_luo_shi = "the Wan-Luo-Shi estimators", x$version_used)
  cat("Imputed mean and SD under ", x$scenario, ", n = ", x$n, ", ", version,
      ", quantile type ", x$type, "\n", sep = "")
  est <- c(mu_hat = x$mu_hat, sigma_hat = x$sigma_hat, se_naive = x$se_naive,
           se_mu = x$se_mu, se_sigma = x$se_sigma)
  print(est, digits = digits, ...)
  cat("Factors: C_mu(n) = ", format(x$c_n, digits = digits),
      ", C_sigma(n) = ", format(x$c_sigma, digits = digits),
      ". The weights are $w_mean and $w_sd; every field with str() or unclass().\n", sep = "")
  invisible(x)
}
