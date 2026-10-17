# ==============================================================================
# File: sdf.R
# Purpose: The standard deviation factors (SDFs) C_mu(n), C_sigma(n) and the inflation
#          factors kappa_mu(n), kappa_sigma(n) of the mean and SD imputed from a
#          five-number summary under normality, from the exact moments of the normal
#          order statistics computed by osmoments::os_moments().
#
# EXPORTED FUNCTIONS:
#   - imputed_moments()   Var/Cov of the imputed mean and SD for a scenario and method
#   - c_mu()              C_mu(n): SD(mu-hat) = C_mu(n) sigma
#   - c_sigma()           C_sigma(n): SD(sigma-hat) = C_sigma(n) sigma
#   - inflation_factors() kappa_mu(n) = C_mu(n) sqrt(n), kappa_sigma(n) = C_sigma(n) / sqrt(1 - b_S(n)^2)
#
# INTERNAL HELPERS:
#   - .check_n_integer()    n is one finite integer within R's integer range
#   - .check_scenario_n()   the smallest sample size of each scenario
#   - .check_exact_rank_n() the exact-rank sample sizes of scenario_ranks()
#   - .b_s()                the bias factor b_S(n) = E[S]/sigma of the sample SD
#   - .var_s()              Var(S)/sigma^2 = 1 - b_S(n)^2, computed without cancellation

#################### EXPORTED MAIN FUNCTIONS ####################

#' Second moments of the imputed mean and SD from a five-number summary
#'
#' Computes \eqn{\mathrm{Var}(\hat\mu)/\sigma^2}, \eqn{\mathrm{Var}(\hat\sigma)/\sigma^2},
#' \eqn{\mathrm{Cov}(\hat\mu,\hat\sigma)/\sigma^2} and \eqn{C(n) = \mathrm{SD}(\hat\mu)/\sigma}
#' for linear estimators \eqn{\hat\mu = w_\mu' Y}, \eqn{\hat\sigma = w_\sigma' Y} of the
#' reported values \eqn{Y}, under normality, using exact order-statistic
#' moments. The identity used is \eqn{\mathrm{Var}(w'Y) = \sigma^2 w'B_Yw}
#' (Balakrishnan et al. 2022, eqs. 49-51), where \eqn{B_Y} is the covariance
#' matrix of the standardized reported values. The default,
#' \code{method_mean = "blue"} with \code{method_sd = "blue"}, is the BLUE, as
#' in every function of the package; \code{method_mean = "luo"} with
#' \code{method_sd = NULL} gives the estimators of Wan, Luo and Shi.
#'
#' @param n Integer sample size: at least 3 under S1 and S2, 5 under S3 and S4,
#'   and 7 under S5, under the default quantile definition. Any such \code{n} is
#'   computed exactly; see Details. Another definition may need a larger \eqn{n}:
#'   where two reported values are the same order statistic (type 3 at
#'   \eqn{n = 5} under S3, for example) the BLUE is not defined and is refused.
#' @param scenario Character. One of \code{"S1"} (min, median, max),
#'   \code{"S2"} (three quartiles), \code{"S3"} (all five). McGrath 2023 labels;
#'   note Wan 2014's C2/C3 are S3/S2 respectively. Also \code{"S4"} (the
#'   quartiles with the 10th and 90th percentiles) and \code{"S5"} (all five with
#'   the two percentiles), the names of Wu, Shi and Tong (2026) for Scenarios 5 and 7 of Balakrishnan
#'   et al. (2022); only the BLUE (\code{"blue"}) is defined for these two.
#' @param method_mean Character, NULL or NA. One of \code{"blue"} (default),
#'   \code{"luo"} or \code{"wan"}, the mean estimator presets of
#'   \code{\link{estimator_weights}}. Ignored if \code{w_mean} is given.
#'   \code{NULL} or \code{NA} skips the mean part.
#' @param method_sd Character, NULL or NA. One of \code{"blue"} (default),
#'   \code{"wan"}, \code{"shi"} (S3 only). \code{NULL} pairs Wan's SD under S1
#'   and S2 with Shi's under S3, the pairing of the Wan-Luo-Shi version (Shi 2020
#'   sec. 2 recommends Wan's SD for S1/S2). Ignored if \code{w_sd} given. Set to
#'   \code{NA} to skip the SD part.
#' @param w_mean,w_sd Optional numeric weight vectors overriding the presets
#'   (length = number of reported values, in increasing order).
#' @param moments Optional precomputed result of \code{osmoments::os_moments(n, ranks)}
#'   for the standard normal, at this \code{n} and at the ranks that the reported
#'   values use (at \code{n = 4Q+1} under type 7 these are
#'   \code{scenario_ranks(n, scenario)}), or a list with its components
#'   \code{alpha} and \code{B} alone. The sample size, the distribution, the
#'   ranks, the location and scale and the shape of \code{alpha} and \code{B} are
#'   checked where the object carries them; a bare list is checked for its size
#'   only. The package keeps the moments it computes itself for the session, so
#'   repeated calls at one sample size integrate once without this argument.
#' @param type The definition of the sample quantiles, an integer from 1 to 9 as
#'   in \code{\link[stats]{quantile}}; 7, R's default, unless a study states
#'   another.
#'
#' @return A list with components:
#'   \code{var_mean} \eqn{= w_\mu'Bw_\mu}, \code{var_sd} \eqn{= w_\sigma'Bw_\sigma},
#'   \code{cov} \eqn{= w_\mu'Bw_\sigma} (all in units of \eqn{\sigma^2});
#'   \code{c_n} \eqn{= \sqrt{w_\mu'Bw_\mu}}; \code{w_mean}, \code{w_sd};
#'   \code{gaps}: named vector of unbiasedness gaps
#'   (\code{mean_sum} \eqn{= w_\mu'1 - 1}, \code{mean_alpha} \eqn{= w_\mu'\alpha}
#'   = bias of \eqn{\hat\mu} in \eqn{\sigma} units,
#'   \code{sd_sum} \eqn{= w_\sigma'1}, \code{sd_alpha} \eqn{= w_\sigma'\alpha - 1}
#'   = relative bias of \eqn{\hat\sigma}); \code{ranks} (the order statistics the
#'   reported values use), \code{n}, \code{scenario}, \code{type}.
#'
#' @details
#' All second moments are proportional to the true \eqn{\sigma^2}; nothing here
#' plugs in an estimate of \eqn{\sigma}. In particular
#' \eqn{\mathrm{SD}(\hat\mu) = C(n)\sigma}, and the plug-in standard error is
#' \eqn{C(n)\hat\sigma}, not \eqn{\hat\sigma/\sqrt{n}}.
#' Each reported value is an order statistic (the minimum, the maximum) or a
#' sample quantile, which under every definition of
#' \code{\link[stats]{quantile}} is a fixed linear combination
#' \eqn{(1-h)X_{(j)} + hX_{(j+1)}} of two adjacent order statistics (Hyndman and
#' Fan 1996). So \eqn{Y = W'X}, with \eqn{W} from
#' \code{osmoments::os_quantile_weights()}, and the means and covariance matrix
#' of the standardized \eqn{Y} are \eqn{\alpha_Y = W'\alpha} and
#' \eqn{B_Y = W'BW}, exact at every \eqn{n}. At \eqn{n = 4Q+1} under type 7 the
#' quartiles are the order statistics at ranks \eqn{Q+1, 2Q+1, 3Q+1}; at other
#' \eqn{n}, or under another definition, they are interpolated, and the result
#' depends on the definition.
#' Assumptions: i.i.d. normal sample; the reported values follow definition
#' \code{type}.
#'
#' @examples
#' im <- imputed_moments(25, "S3")            # the BLUE of the mean and of the SD, the default
#' im$c_n                                     # SD(mu-hat)/sigma
#' imputed_moments(25, "S3", "luo", NULL)$c_n # Luo's mean, with Shi's SD
#' sqrt(1/25)                                 # naive sigma/sqrt(n) factor, for contrast
#'
#' @references
#' Balakrishnan N, Rychtar J, Taylor D, Walter SD (2022). Unified approach to
#' optimal estimation of mean and standard deviation from sample summaries.
#' Statistical Methods in Medical Research 31(11), 2087-2103. (eqs. 49-51)
#' Luo D, Wan X, Liu J, Tong T (2018). Optimally estimating the sample mean from
#' the sample size, median, mid-range, and/or mid-quartile range. Statistical
#' Methods in Medical Research 27(6), 1785-1805.
#' Shi J, Luo D, Weng H, Zeng XT, Lin L, Chu H, Tong T (2020). Optimally
#' estimating the sample standard deviation from the five-number summary.
#' Research Synthesis Methods 11(5), 641-654.
#' Wan X, Wang W, Liu J, Tong T (2014). Estimating the sample mean and standard
#' deviation from the sample size, median, range and/or interquartile range.
#' BMC Medical Research Methodology 14, 135.
#' Yang X, Hutson AD, Wang D (2022). A generalized BLUE approach for combining
#' location and scale information in a meta-analysis. Journal of Applied
#' Statistics 49(15), 3846-3867. (eq. 22)
#' Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard
#' deviation from a five-number summary. Manuscript.
#'
#' @export
imputed_moments <- function(n, scenario = c("S3", "S1", "S2", "S4", "S5"),
                            method_mean = "blue", method_sd = "blue",
                            w_mean = NULL, w_sd = NULL, moments = NULL, type = 7) {
  scenario <- match.arg(scenario)
  .check_scenario_n(n, scenario)
  # the means and covariance matrix of the standardized reported values, Y = W'X
  rm <- .reported_moments(n, scenario, type, moments)
  moments <- rm$moments; ranks <- rm$ranks
  alpha <- rm$alpha; B <- rm$B

  # method_sd = NULL: the Wan-Luo-Shi pairing, Shi for S3, Wan for S1/S2 (see @param method_sd)
  if (is.null(method_sd)) method_sd <- if (scenario == "S3") "shi" else "wan"

  if (is.null(w_mean) && !is.null(method_mean) && !is.na(method_mean))
    w_mean <- estimator_weights(n, scenario, "mean", method_mean, moments = moments, type = type)
  if (is.null(w_sd) && !is.null(method_sd) && !is.na(method_sd))
    w_sd <- estimator_weights(n, scenario, "sd", method_sd, moments = moments, type = type)

  qf <- function(u, v) as.numeric(t(u) %*% B %*% v)  # u' B_Y v, sigma^2 units
  var_mean <- if (!is.null(w_mean)) qf(w_mean, w_mean) else NA_real_
  var_sd   <- if (!is.null(w_sd))   qf(w_sd, w_sd)     else NA_real_
  cov_ms   <- if (!is.null(w_mean) && !is.null(w_sd)) qf(w_mean, w_sd) else NA_real_

  gaps <- c(
    mean_sum   = if (!is.null(w_mean)) sum(w_mean) - 1        else NA_real_,
    mean_alpha = if (!is.null(w_mean)) sum(w_mean * alpha)     else NA_real_,
    sd_sum     = if (!is.null(w_sd))   sum(w_sd)               else NA_real_,
    sd_alpha   = if (!is.null(w_sd))   sum(w_sd * alpha) - 1   else NA_real_
  )

  list(n = n, scenario = scenario, ranks = ranks, type = type,
       w_mean = w_mean, w_sd = w_sd,
       var_mean = var_mean, var_sd = var_sd, cov = cov_ms,
       c_n = sqrt(var_mean), gaps = gaps)
}

#' The scalar C_mu(n): standard deviation of the imputed mean in units of sigma
#'
#' \eqn{C_\mu(n) = \sqrt{w_\mu' B w_\mu}} so that \eqn{\mathrm{SD}(\hat\mu) = C_\mu(n)\,\sigma}
#' (the true standard deviation of the estimator, a function of the unknown
#' population \eqn{\sigma}). The plug-in standard error is
#' \eqn{\mathrm{SE}(\hat\mu) = C_\mu(n)\,\hat\sigma}; its reference distribution is
#' the subject of a separate study. Exact (numerical integration); for the
#' closed-form shortcut formula see \code{\link{c_mu_shortcut}}.
#'
#' @inheritParams imputed_moments
#' @param method Mean estimator preset; see \code{\link{estimator_weights}}. The
#'   default \code{"blue"} is the BLUE, the default estimator version of Wu, Shi and
#'   Tong (2026);
#'   \code{"luo"} gives the factor of Luo's shortcut mean estimator.
#' @return Numeric scalar \eqn{C(n)}.
#' @examples
#' c_mu(25, "S2")                 # exact C(n) for the BLUE of the mean under S2
#' c_mu(25, "S2", "luo")          # the same for Luo's estimator
#' c_mu(25, "S2") * sqrt(25)      # inflation factor over sigma/sqrt(n)
#' @references
#' Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard
#' deviation from a five-number summary. Manuscript.
#' @export
c_mu <- function(n, scenario = c("S3", "S1", "S2", "S4", "S5"), method = "blue", moments = NULL,
                 type = 7) {
  scenario <- match.arg(scenario)
  im <- imputed_moments(n, scenario, method_mean = method, method_sd = NA,
                        moments = moments, type = type)
  im$c_n
}

#' The scalar C_sigma(n): standard deviation of the imputed SD in units of sigma
#'
#' \eqn{C_\sigma(n) = \sqrt{w_\sigma' B w_\sigma}} so that
#' \eqn{\mathrm{SD}(\hat\sigma) = C_\sigma(n)\,\sigma}, the true standard
#' deviation of the imputed SD (a function of the unknown population
#' \eqn{\sigma}); the plug-in standard error is
#' \eqn{\mathrm{SE}(\hat\sigma) = C_\sigma(n)\,\hat\sigma}. The factor is exact
#' (numerical integration of the order-statistic moments), and
#' \code{\link{c_sigma_shortcut}} gives its approximation formula. This is the
#' SD companion of \code{\link{c_mu}}.
#'
#' @inheritParams imputed_moments
#' @param method SD estimator preset (see \code{\link{estimator_weights}}). The
#'   default \code{"blue"} is the BLUE of \eqn{\sigma}, the default estimator
#'   version of Wu, Shi and Tong (2026); \code{NULL} selects the SD estimator of Wan
#'   (S1 and S2) or Shi (S3).
#' @return Numeric scalar \eqn{C_\sigma(n)}.
#' @details Assumptions as in \code{\link{imputed_moments}}: iid normal data, the
#'   reported values following definition \code{type}. The value is
#'   \code{sqrt(imputed_moments(...)$var_sd)}.
#' @examples
#' c_sigma(25, "S2")                    # exact C_sigma(n) for the BLUE of sigma under S2
#' c_sigma(25, "S2", method = NULL)     # the same for Wan's estimator
#' c_sigma(25, "S2") * sqrt(2 * 25)     # inflation over the large-sample sigma/sqrt(2n)
#' @references
#' Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard
#' deviation from a five-number summary. Manuscript.
#' @export
c_sigma <- function(n, scenario = c("S3", "S1", "S2", "S4", "S5"), method = "blue",
                    moments = NULL, type = 7) {
  scenario <- match.arg(scenario)
  im <- imputed_moments(n, scenario, method_mean = NA, method_sd = method,
                        moments = moments, type = type)
  sqrt(im$var_sd)
}

#' Inflation factors of the imputed mean and SD
#'
#' The inflation factors (IF) of Wu, Shi and Tong (2026): the ratio of the
#' standard deviation of the imputed mean to that of the sample mean,
#' \eqn{\kappa_\mu(n) = \mathrm{SD}(\hat\mu)/\mathrm{SD}(\bar X) = C_\mu(n)\sqrt n},
#' and of the imputed SD to that of the sample SD \eqn{S},
#' \eqn{\kappa_\sigma(n) = \mathrm{SD}(\hat\sigma)/\mathrm{SD}(S) = C_\sigma(n)/\sqrt{1 - b_S(n)^2}},
#' where \eqn{b_S(n) = \sqrt{2/(n-1)}\,\Gamma(n/2)/\Gamma\{(n-1)/2\}} is the bias
#' factor of \eqn{S}, \eqn{E[S] = b_S(n)\sigma}. The factors \eqn{C_\mu(n)} and
#' \eqn{C_\sigma(n)} are the exact ones of \code{\link{c_mu}} and
#' \code{\link{c_sigma}}, from one exact moment set.
#'
#' @inheritParams imputed_moments
#' @param version \code{"blue"} (default) for the BLUE of both, or
#'   \code{"wan_luo_shi"} for the estimators of Wan, Luo and Shi (Luo's mean; Wan's
#'   SD under S1 and S2, Shi's under S3). Under S4 and S5 only \code{"blue"} is
#'   defined.
#' @return Named numeric vector: \code{kappa_mu}, \code{kappa_sigma}.
#' @details Assumptions as in \code{\link{imputed_moments}}: iid normal data, the
#'   reported values following definition \code{type}. An IF of 1.2 means a 20\%
#'   larger SD than that of the full-sample statistic.
#' @examples
#' inflation_factors(41, "S2")                       # the BLUE: about 1.076 and 1.651
#' inflation_factors(41, "S2", version = "wan_luo_shi")
#' inflation_factors(40, "S2")                       # an even n: interpolated quartiles, type 7
#' inflation_factors(21, "S4")                       # quartiles with the 10th and 90th percentiles
#' @references
#' Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard
#' deviation from a five-number summary. Manuscript.
#' @export
inflation_factors <- function(n, scenario = c("S3", "S1", "S2", "S4", "S5"),
                              version = c("blue", "wan_luo_shi"), type = 7) {
  scenario <- match.arg(scenario); version <- match.arg(version)
  .check_scenario_n(n, scenario)
  if (version == "wan_luo_shi" && scenario %in% c("S4", "S5"))
    stop("the Wan-Luo-Shi estimators are not defined for ", scenario, "; use version = \"blue\"")
  # one exact moment set serves both factors
  mo <- .reported_moments(n, scenario, type)$moments
  m_mean <- if (version == "blue") "blue" else "luo"
  m_sd   <- if (version == "blue") "blue" else if (scenario == "S3") "shi" else "wan"
  cm <- c_mu(n, scenario, method = m_mean, moments = mo, type = type)
  cs <- c_sigma(n, scenario, method = m_sd, moments = mo, type = type)
  # SD(xbar) = sigma / sqrt(n) and, under normality, SD(S) = sigma sqrt(1 - b_S(n)^2)
  c(kappa_mu = cm * sqrt(n), kappa_sigma = cs / sqrt(.var_s(n)))
}

#################### INTERNAL HELPER FUNCTIONS ####################

#' n is one finite integer within R's integer range
#'
#' Checked before any as.integer(), which would turn Inf or a value above
#' .Machine$integer.max into NA with a warning and a later obscure error.
#' @keywords internal
#' @noRd
.check_n_integer <- function(n) {
  if (length(n) != 1 || !is.numeric(n) || !is.finite(n) || n != round(n) || n > .Machine$integer.max)
    stop("n must be a single finite integer (at most ", .Machine$integer.max, "); got: ",
         paste(n, collapse = ", "))
  invisible(TRUE)
}

#' Validate the (n, scenario) combination
#' @keywords internal
#' @noRd
.check_scenario_n <- function(n, scenario) {
  .check_n_integer(n)
  # the smallest n at which the reported values are distinct linear combinations
  # of the order statistics (under type 7): as many order statistics as values
  nmin <- c(S1 = 3, S2 = 3, S3 = 5, S4 = 5, S5 = 7)[[scenario]]
  if (n < nmin)
    stop(scenario, " requires n >= ", nmin, "; got n = ", n)
  invisible(TRUE)
}

#' The exact-rank sample sizes of scenario_ranks()
#' @keywords internal
#' @noRd
.check_exact_rank_n <- function(n, scenario) {
  .check_n_integer(n)
  if (scenario == "S1") {
    if (n < 3 || n %% 2 == 0)
      stop("S1 reports a single median rank at odd n >= 3; got n = ", n)
  } else if (scenario %in% c("S4", "S5")) {
    if (n < 21 || (n - 1) %% 20 != 0)
      stop(scenario, " reports single ranks at n = 20m+1 with integer m >= 1; got n = ", n)
  } else {
    if (n < 5 || (n - 1) %% 4 != 0)
      stop(scenario, " reports single ranks at n = 4Q+1 with integer Q >= 1; got n = ", n)
  }
  invisible(TRUE)
}

#' Bias factor of the sample SD under normality, and the variance of S
#'
#' \eqn{b_S(n) = E[S]/\sigma = \sqrt{2/(n-1)}\,\Gamma(n/2)/\Gamma\{(n-1)/2\}
#' = \sqrt{2\pi/(n-1)}/B\{(n-1)/2, 1/2\}}, through \code{lbeta()}, which keeps
#' its digits at large \eqn{n}, where the difference of two \code{lgamma()} values
#' gave \eqn{b_S > 1} at \eqn{n = 10^8}. \code{.var_s(n)}
#' \eqn{= \mathrm{Var}(S)/\sigma^2 = 1 - b_S(n)^2}, computed as
#' \eqn{-\mathrm{expm1}\{\log b_S(n)^2\}} so that it stays positive when
#' \eqn{b_S(n)} is within rounding of 1; it is \eqn{1/(2n)} to first order.
#' @keywords internal
#' @noRd
.b_s <- function(n) sqrt(2 * pi / (n - 1)) / exp(lbeta((n - 1) / 2, 0.5))
.var_s <- function(n) -expm1(log(2 * pi / (n - 1)) - 2 * lbeta((n - 1) / 2, 0.5))
