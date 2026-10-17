# ==============================================================================
# File: se_imputed_df.R
# Purpose: The data-frame interface for a meta-analysis: one row per study, the
#          reporting scenario read off the columns present, the imputed mean and
#          SD with their naive and exact plug-in standard errors appended; and the
#          two-group effect per row (the difference of means or the log
#          variability ratio) with its naive and exact standard errors, ready for
#          pool_imputed() and for metafor::rma().
#
# EXPORTED FUNCTIONS:
#   - se_imputed_df()    every study of a data frame through se_imputed()
#   - se_imputed_diff()  two groups per row: the effect (measure "MD" or "lnVR") and its SEs
#
# INTERNAL HELPERS:
#   - .read_scenario()   which scenario a row's non-missing values give
#   - .range_sd()           the SD imputed from the range alone, with its exact factor

#################### EXPORTED MAIN FUNCTIONS ####################

#' Imputed means and SDs with their standard errors for every study of a data frame
#'
#' The data-frame interface of \code{\link{se_imputed}} for a meta-analysis: one
#' row per study, with the columns that such a data set usually carries (the
#' sample size, the minimum, the first quartile, the median, the third quartile,
#' the maximum, and the sample mean and SD where a study gives them). For each
#' row the function reads which values are present, decides the data-reporting
#' scenario, imputes the mean and the SD where they are not given, and appends
#' them with their standard errors: the naive \eqn{\hat\sigma/\sqrt n} that
#' treats the imputed mean as a sample mean, and the plug-in standard errors
#' \eqn{C_\mu(n)\hat\sigma} and \eqn{C_\sigma(n)\hat\sigma} built from the
#' factors of Wu, Shi and Tong (2026).
#'
#' @param data A data frame, one row per study.
#' @param n,min,q1,med,q3,max,mean,sd The names of the columns holding the sample
#'   size, the minimum, the first quartile, the median, the third quartile, the
#'   maximum, and the sample mean and SD, each with \code{suffix} appended. An
#'   argument set to \code{NULL}, or naming a column that \code{data} does not
#'   have, counts as absent for every row. A missing value (\code{NA}) counts as
#'   absent for that row. The columns must be numeric; a column that is all
#'   \code{NA} may be logical, as \code{read.csv()} reads such a column.
#' @param suffix A string appended to every column name above and to every output
#'   column, for data sets that carry two groups side by side (\code{".g1"},
#'   \code{".g2"}). The default is the empty string.
#' @param version \code{"blue"} (default) or \code{"wan_luo_shi"}, as in
#'   \code{\link{se_imputed}}.
#' @param type The definition of the reported sample quantiles, 1 to 9 as in
#'   \code{\link[stats]{quantile}}; 7, R's default. One value for every row.
#' @param prefer What to use for a row that has both a summary and a sample
#'   mean and SD: \code{"sample"} (default), the sample statistics, or
#'   \code{"imputed"}.
#' @param append \code{TRUE} (default) returns \code{data} with the output
#'   columns appended; \code{FALSE} returns the output columns alone.
#' @return A data frame. The output columns, each with \code{suffix} appended, are
#'   \code{scenario} (\code{"S1"}, \code{"S2"}, \code{"S3"}, \code{"sample"} or
#'   \code{"mean_range"}, see Details), \code{mu_hat}, \code{sigma_hat},
#'   \code{se_naive} (\eqn{\hat\sigma/\sqrt n}), \code{se_mu}
#'   (\eqn{C_\mu(n)\hat\sigma} for an imputed mean, \eqn{\hat\sigma/\sqrt n} for a
#'   sample mean), \code{se_sigma} (\eqn{C_\sigma(n)\hat\sigma} for an imputed SD,
#'   and for a sample SD \eqn{\mathrm{sd}\sqrt{1 - b_S(n)^2}}, the standard
#'   deviation of a sample SD under normality with the SD plugged in),
#'   \code{version_used} (\code{"blue"}, \code{"wan_luo_shi"}, or \code{"sample"}
#'   when nothing is imputed) and \code{exact} (\code{TRUE} for a row with an
#'   imputed value, whose factors are computed exactly at its \eqn{n}, \code{NA}
#'   for a \code{"sample"} row).
#' @details
#' The scenario of a row is read off its non-missing values, in this order. A
#' sample mean and SD give \code{"sample"}, the two used as they are, unless
#' \code{prefer = "imputed"} and the row also has a complete summary. All five
#' values give \code{"S3"}; the first quartile, the median and the third
#' quartile give \code{"S2"} (a single extreme beside them is ignored); the
#' minimum, the median and the maximum give \code{"S1"}. A mean with the
#' minimum and the maximum, and no SD, gives \code{"mean_range"} when no
#' summary is complete: the mean is kept as the sample mean, and the SD is
#' imputed from the two extremes alone, as the range divided by its expected
#' value in a normal sample of size \eqn{n} (the BLUE of \eqn{\sigma} from the
#' extremes, under \code{version = "blue"}) or by the closed-form divisor
#' \eqn{\xi(n)} of the Wan-Luo-Shi version (under \code{"wan_luo_shi"}), with
#' \code{se_sigma} \eqn{= C_\sigma(n)\hat\sigma} from the exact moments of the
#' two extremes. These are the SD estimators of \code{"S1"} for a sample whose
#' median sits at the center of the order statistics, where the median carries
#' no weight, computed here without a median, so the quantile definition
#' \code{type} plays no part. Such a row has \code{se_naive} and \code{se_mu}
#' both equal to \eqn{\hat\sigma/\sqrt n}, the standard error of a sample mean
#' with the imputed SD plugged in, and needs \eqn{n \ge 3}. A mean outside its
#' range gets \code{NA}, as an error in the row. A mean without an SD beside a
#' complete summary is ignored, and the summary gives both estimates.
#'
#' A row with none of these, with a missing, infinite or non-integer sample
#' size, with a non-finite mean, with a non-finite, negative or zero SD, with a
#' summary whose values are all equal or a range of zero (its imputed SD is
#' zero), with values out of order or a mean outside its reported range, gets
#' \code{NA} in every output column and
#' is named in one warning with its reason; the function does not stop, so that
#' the whole table comes back. A standard error of zero would give the study an
#' infinite inverse-variance weight, which is why such a row is set to \code{NA}
#' rather than kept. Output columns that an earlier call left in \code{data} are
#' replaced, so a data frame can be run again after a value, the \code{version}
#' or the \code{type} changes.
#'
#' Every imputed row is computed exactly at its own \eqn{n}, as in
#' \code{\link{se_imputed}}: at a sample size not of the form \eqn{4Q + 1} the
#' reported quartiles are interpolated under definition \code{type}.
#'
#' Assumptions: independent normal data within each study. The standard errors
#' are plug-in quantities, \eqn{\sigma} replaced by \eqn{\hat\sigma}; their
#' reference distribution is not Student's t and is not supplied by this version.
#' A study whose data are visibly skewed (the quartiles far from symmetric around
#' the median) should be screened before the normal-model estimators are applied,
#' as the source papers of those estimators recommend.
#' @examples
#' studies <- data.frame(n = c(25, 30, 41, 25, 35),
#'                       min = c(3, NA, 1, NA, 2), q1 = c(10, 10, 9, NA, NA),
#'                       med = c(12, 12, 12, NA, NA), q3 = c(15, 15, 16, NA, NA),
#'                       max = c(30, NA, 30, NA, 28),
#'                       mean = c(NA, NA, NA, 13, 12), sd = c(NA, NA, NA, 3, NA))
#' se_imputed_df(studies)   # S3, S2, S1, a sample mean and SD, a mean with the range
#' se_imputed_df(studies, version = "wan_luo_shi", append = FALSE)
#' # columns named as a spreadsheet often names them
#' names(studies) <- c("N", "Min", "Q1", "Median", "Q3", "Max", "Mean", "SD")
#' se_imputed_df(studies, n = "N", min = "Min", q1 = "Q1", med = "Median", q3 = "Q3",
#'               max = "Max", mean = "Mean", sd = "SD", append = FALSE)
#' @seealso \code{\link{se_imputed}} for one study, \code{\link{se_imputed_diff}}
#'   for two groups per row, \code{\link{pool_imputed}} to pool the studies.
#' @references
#' Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard
#' deviation from a five-number summary. Manuscript.
#' @export
se_imputed_df <- function(data, n = "n", min = "min", q1 = "q1", med = "med", q3 = "q3",
                          max = "max", mean = "mean", sd = "sd", suffix = "",
                          version = c("blue", "wan_luo_shi"), prefer = c("sample", "imputed"),
                          append = TRUE, type = 7) {
  version <- match.arg(version); prefer <- match.arg(prefer)
  if (!is.data.frame(data)) stop("data must be a data frame, one row per study")
  if (!is.character(suffix) || length(suffix) != 1 || is.na(suffix))
    stop("suffix must be one string")
  # a named column's values, or NULL when the argument is NULL or the column is absent
  col <- function(name) {
    if (is.null(name)) return(NULL)
    if (!is.character(name) || length(name) != 1) stop("a column argument must be one name")
    nm <- paste0(name, suffix)
    if (!nm %in% names(data)) return(NULL)
    v <- data[[nm]]
    if (is.logical(v) && all(is.na(v))) v <- as.numeric(v)   # an all-NA column, read in as logical
    if (!is.numeric(v)) stop("column ", nm, " must be numeric")
    v
  }
  v_n <- col(n); v_min <- col(min); v_q1 <- col(q1); v_med <- col(med)
  v_q3 <- col(q3); v_max <- col(max); v_mean <- col(mean); v_sd <- col(sd)
  K <- nrow(data)
  pick <- function(v, i) if (is.null(v)) NA_real_ else as.numeric(v[i])
  na_r <- rep(NA_real_, K); na_c <- rep(NA_character_, K)   # every column of length K, K = 0 included
  out <- data.frame(scenario = na_c, mu_hat = na_r, sigma_hat = na_r,
                    se_naive = na_r, se_mu = na_r, se_sigma = na_r,
                    version_used = na_c, exact = rep(NA, K), stringsAsFactors = FALSE)
  bad <- character(0)
  for (i in seq_len(K)) {
    ni <- pick(v_n, i)
    vals <- c(min = pick(v_min, i), q1 = pick(v_q1, i), med = pick(v_med, i),
              q3 = pick(v_q3, i), max = pick(v_max, i))
    sample <- c(mean = pick(v_mean, i), sd = pick(v_sd, i))
    sc <- .read_scenario(vals, sample, prefer)
    if (is.na(sc)) { bad <- c(bad, paste0(i, " (no usable summary)")); next }
    if (!is.finite(ni) || ni < 2 || ni != round(ni) || ni > .Machine$integer.max) {
      bad <- c(bad, paste0(i, " (sample size)")); next
    }
    ni <- as.integer(round(ni))
    if (sc == "sample") {
      s <- sample[["sd"]]
      if (!is.finite(sample[["mean"]]) || !is.finite(s)) { bad <- c(bad, paste0(i, " (non-finite mean or SD)")); next }
      if (s < 0) { bad <- c(bad, paste0(i, " (negative SD)")); next }
      if (s == 0) { bad <- c(bad, paste0(i, " (zero SD)")); next }   # a zero SE cannot enter the weights
      # a sample mean and SD: the sample mean's SE, and the SD of the sample SD
      # under normality, sigma sqrt(1 - b_S(n)^2), with the SD plugged in
      out$scenario[i] <- "sample"; out$mu_hat[i] <- sample[["mean"]]; out$sigma_hat[i] <- s
      out$se_naive[i] <- s / sqrt(ni); out$se_mu[i] <- s / sqrt(ni)
      out$se_sigma[i] <- s * sqrt(.var_s(ni))
      out$version_used[i] <- "sample"; out$exact[i] <- NA
      next
    }
    if (sc == "mean_range") {
      # a sample mean with the range and no SD: the mean is kept, and the SD comes from
      # the two extremes alone (.range_sd()), with no median and so no dependence on `type`
      m <- sample[["mean"]]; a <- vals[["min"]]; b <- vals[["max"]]
      if (!is.finite(m)) { bad <- c(bad, paste0(i, " (non-finite mean)")); next }
      if (!is.finite(a) || !is.finite(b)) {
        bad <- c(bad, paste0(i, " (non-finite extremes)")); next
      }
      if (a > b) { bad <- c(bad, paste0(i, " (minimum above maximum)")); next }
      if (a == b) { bad <- c(bad, paste0(i, " (zero range)")); next }   # an imputed SD of zero
      if (m < a || m > b) { bad <- c(bad, paste0(i, " (mean outside the range)")); next }
      if (ni < 3) { bad <- c(bad, paste0(i, " (a mean with the range needs n >= 3)")); next }
      r <- .range_sd(a, b, ni, version)
      out$scenario[i] <- "mean_range"; out$mu_hat[i] <- m; out$sigma_hat[i] <- r$sigma_hat
      # the standard error of a sample mean, with the imputed SD plugged in, under both rules
      out$se_naive[i] <- r$sigma_hat / sqrt(ni); out$se_mu[i] <- r$sigma_hat / sqrt(ni)
      out$se_sigma[i] <- r$c_sigma * r$sigma_hat
      out$version_used[i] <- version; out$exact[i] <- TRUE
      next
    }
    q <- switch(sc,
      S1 = vals[c("min", "med", "max")],
      S2 = vals[c("q1", "med", "q3")],
      S3 = vals[c("min", "q1", "med", "q3", "max")])
    # a constant summary has an imputed SD of zero, so zero standard errors: NA, as a zero sample SD
    if (isTRUE(max(q) == min(q))) { bad <- c(bad, paste0(i, " (constant summary)")); next }
    r <- tryCatch(se_imputed(q, ni, sc, version = version, type = type), error = function(e) e)
    if (inherits(r, "error")) { bad <- c(bad, paste0(i, " (", conditionMessage(r), ")")); next }
    out$scenario[i] <- sc; out$mu_hat[i] <- r$mu_hat; out$sigma_hat[i] <- r$sigma_hat
    out$se_naive[i] <- r$se_naive; out$se_mu[i] <- r$se_mu; out$se_sigma[i] <- r$se_sigma
    out$version_used[i] <- r$version_used; out$exact[i] <- r$exact
  }
  if (length(bad))
    warning("rows with NA output, ", length(bad), " of ", K, ": ", paste(bad, collapse = "; "))
  names(out) <- paste0(names(out), suffix)
  if (!append) return(out)
  # the output columns of an earlier call are replaced, not duplicated: with duplicate names
  # se_imputed_diff() would read the stale ones
  cbind(data[setdiff(names(data), names(out))], out)
}

#' The effect of two groups per study, with its naive and exact standard errors
#'
#' For a data set with two groups side by side (the columns of one group carry
#' one suffix, those of the other group another), runs
#' \code{\link{se_imputed_df}} on each group and appends the effect of the
#' study, group 1 against group 2, with its standard error under the naive and
#' under the exact standard errors: the difference of means
#' (\code{measure = "MD"}, the default) or the log variability ratio
#' (\code{measure = "lnVR"}), in the column names that \code{metafor::rma()}
#' takes as they are, \code{yi} and \code{sei}. \code{\link{pool_imputed}} pools
#' the result in one call.
#'
#' @inheritParams se_imputed_df
#' @param suffixes The two suffixes of the groups' columns, group 1 first. The
#'   default is the layout of the data sets shipped with the metamedian package,
#'   \code{".g1"} and \code{".g2"}.
#' @param ... The column-name arguments of \code{\link{se_imputed_df}} (\code{n},
#'   \code{min}, \code{q1}, \code{med}, \code{q3}, \code{max}, \code{mean},
#'   \code{sd}), passed to both groups.
#' @param measure The effect per study: \code{"MD"} (default), the difference of
#'   the two means, group 1 minus group 2; or \code{"lnVR"}, the log variability
#'   ratio of Nakagawa et al. (2015), the log of the SD of group 1 over the SD
#'   of group 2. See Details.
#' @return \code{data} with the output columns of \code{\link{se_imputed_df}} for
#'   each group (each with its suffix) and three more: \code{yi}, the effect;
#'   \code{sei_naive}, its standard error from the naive standard errors; and
#'   \code{sei}, its standard error from the exact ones. The result carries the
#'   attribute \code{"measure"}, which \code{\link{pool_imputed}} reads for its
#'   heading. A row that \code{se_imputed_df()} could not use in either group has
#'   \code{NA} in the three columns.
#' @details
#' The two groups are independent samples, so the variances add.
#'
#' Under \code{"MD"}, \eqn{y_i = \hat\mu_{i1} - \hat\mu_{i2}},
#' \code{sei_naive} \eqn{= (\mathrm{se\_naive}_{i1}^2 + \mathrm{se\_naive}_{i2}^2)^{1/2}}
#' and \code{sei} \eqn{= (\mathrm{se\_mu}_{i1}^2 + \mathrm{se\_mu}_{i2}^2)^{1/2}}.
#'
#' Under \code{"lnVR"}, \eqn{y_i = \log(\hat\sigma_{i1}/\hat\sigma_{i2})}, and the
#' variance of \eqn{y_i} is approximated by the sum of one term per group, to
#' first order in each group's SD (the delta method,
#' \eqn{\mathrm{Var}(\log X) \approx \mathrm{Var}(X)/(\mathrm{E}X)^2}). The
#' naive variance gives every group the large-sample term of the log of a sample
#' SD of \eqn{n} normal observations, \eqn{1/\{2(n - 1)\}}, so
#' \code{sei_naive} \eqn{= \{1/(2(n_{i1} - 1)) + 1/(2(n_{i2} - 1))\}^{1/2}}. The
#' exact variance keeps that term for a group that gives its sample SD and uses
#' \eqn{C^2_\sigma(n)}, the exact variance of the imputed SD in units of
#' \eqn{\sigma^2}, for a group whose SD is imputed (from the quartiles, the
#' extremes, or the range alone), which the group's columns carry as
#' \code{(se_sigma / sigma_hat)^2}; \code{sei} is the square root of the sum of
#' the two terms. These are the two variances of Wu, Shi and Tong (2026): exact
#' in the factor \eqn{C_\sigma(n)}, first-order in the log, so neither is the
#' finite-sample variance of the log ratio, and the intervals of
#' \code{\link{pool_imputed}} are normal approximations. The small-sample bias
#' correction of Nakagawa et al. (2015), which
#' \code{metafor::escalc(measure = "VR")} adds to \eqn{y_i}, is not applied, as
#' in Wu, Shi and Tong (2026). Both analyses use the plain log ratio.
#' @examples
#' two <- data.frame(n.g1 = c(25, 30), q1.g1 = c(10, 11), med.g1 = c(12, 13), q3.g1 = c(15, 16),
#'                   n.g2 = c(29, 40), q1.g2 = c(8, 9), med.g2 = c(10, 11), q3.g2 = c(13, 14))
#' d <- se_imputed_diff(two)
#' d[, c("yi", "sei_naive", "sei")]
#' pool_imputed(d)                                    # the pooled difference, both SEs
#' v <- se_imputed_diff(two, measure = "lnVR")       # the log ratio of the SDs
#' v[, c("yi", "sei_naive", "sei")]
#' # columns named otherwise, two groups with the suffixes "_t" and "_c"
#' spread <- data.frame(N_t = 25, Min_t = 3, Median_t = 12, Max_t = 30,
#'                      N_c = 29, Min_c = 2, Median_c = 10, Max_c = 25)
#' se_imputed_diff(spread, suffixes = c("_t", "_c"), n = "N", min = "Min",
#'                 med = "Median", max = "Max")[, c("yi", "sei_naive", "sei")]
#' \dontrun{
#' metafor::rma(yi, sei = sei, data = d)          # the pooled difference with the exact SEs
#' metafor::rma(yi, sei = sei_naive, data = d)    # the same with the naive SEs, for comparison
#' }
#' @seealso \code{\link{se_imputed_df}}, \code{\link{se_imputed}},
#'   \code{\link{pool_imputed}}.
#' @references
#' Nakagawa S, Poulin R, Mengersen K, Reinhold K, Engqvist L, Lagisz M, Senior AM
#' (2015). Meta-analysis of variation: ecological and evolutionary applications
#' and beyond. Methods in Ecology and Evolution 6(2), 143-152.
#'
#' Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard
#' deviation from a five-number summary. Manuscript.
#' @export
se_imputed_diff <- function(data, suffixes = c(".g1", ".g2"), ..., measure = c("MD", "lnVR"),
                            version = c("blue", "wan_luo_shi"), prefer = c("sample", "imputed"),
                            type = 7) {
  version <- match.arg(version); prefer <- match.arg(prefer); measure <- match.arg(measure)
  if (!is.character(suffixes) || length(suffixes) != 2 || anyNA(suffixes) ||
      suffixes[1] == suffixes[2])
    stop("suffixes must be two different strings, group 1 first")
  d <- se_imputed_df(data, ..., suffix = suffixes[1], version = version, prefer = prefer,
                     append = TRUE, type = type)
  d <- se_imputed_df(d, ..., suffix = suffixes[2], version = version, prefer = prefer,
                     append = TRUE, type = type)
  g <- function(stem, s) d[[paste0(stem, s)]]
  if (measure == "MD") {
    d$yi        <- g("mu_hat", suffixes[1]) - g("mu_hat", suffixes[2])
    d$sei_naive <- sqrt(g("se_naive", suffixes[1])^2 + g("se_naive", suffixes[2])^2)
    d$sei       <- sqrt(g("se_mu", suffixes[1])^2 + g("se_mu", suffixes[2])^2)
  } else {
    # the sample-size column: the argument n of se_imputed_df() as the dots give it, by
    # name, by position or by partial match, matched exactly as se_imputed_df() matched it
    dots <- list(...)
    m <- match.call(se_imputed_df, as.call(c(list(as.name("se_imputed_df"), quote(data)), dots)))
    n_name <- if ("n" %in% names(m)) m[["n"]] else "n"
    n1 <- if (is.null(n_name)) NULL else g(n_name, suffixes[1])
    n2 <- if (is.null(n_name)) NULL else g(n_name, suffixes[2])
    if (is.null(n1) || is.null(n2))
      stop("measure = \"lnVR\" needs the sample size of each group (the column n)")
    # the rows usable in both groups; the others get NA in the three columns, as under "MD"
    ok <- !is.na(g("scenario", suffixes[1])) & !is.na(g("scenario", suffixes[2]))
    # one group's term of the variance of the log SD: 1 / (2 (n - 1)) for a sample
    # SD, C_sigma(n)^2 = (se_sigma / sigma_hat)^2 for an imputed one
    term <- function(s) {
      ifelse(g("scenario", s) == "sample", 1 / (2 * (as.numeric(g(n_name, s)) - 1)),
             (g("se_sigma", s) / g("sigma_hat", s))^2)
    }
    yi <- sn <- se <- rep(NA_real_, nrow(d))
    # the difference of the logs, finite where the ratio itself would overflow
    yi[ok] <- log(g("sigma_hat", suffixes[1])[ok]) - log(g("sigma_hat", suffixes[2])[ok])
    sn[ok] <- sqrt(1 / (2 * (as.numeric(n1)[ok] - 1)) + 1 / (2 * (as.numeric(n2)[ok] - 1)))
    se[ok] <- sqrt(term(suffixes[1])[ok] + term(suffixes[2])[ok])
    d$yi <- yi; d$sei_naive <- sn; d$sei <- se
  }
  attr(d, "measure") <- measure
  d
}

#################### INTERNAL HELPER FUNCTIONS ####################

#' Which scenario a row's non-missing values give
#'
#' \code{vals} is the named vector (min, q1, med, q3, max) with NA for an absent
#' value, \code{sample} the named vector (mean, sd), \code{prefer} "sample" or
#' "imputed". Returns "S1", "S2", "S3", "sample", "mean_range" or NA, in the
#' order the help of se_imputed_df() states.
#' @keywords internal
#' @noRd
.read_scenario <- function(vals, sample, prefer) {
  has <- !is.na(vals)
  got <- !is.na(sample)
  sc <- if (all(has)) "S3" else
        if (all(has[c("q1", "med", "q3")])) "S2" else
        if (all(has[c("min", "med", "max")])) "S1" else NA_character_
  if (all(got) && (prefer == "sample" || is.na(sc))) return("sample")
  # a mean with the minimum and the maximum, no SD, and no complete summary
  if (is.na(sc) && got[["mean"]] && !got[["sd"]] && has[["min"]] && has[["max"]]) return("mean_range")
  sc
}

#' The SD imputed from the range alone, with its exact standard deviation factor
#'
#' For a row that gives a mean with the minimum and the maximum and no SD. The
#' estimator is the range over a divisor: under "blue" the expected range of a
#' standard normal sample of size n, alpha_(n) - alpha_(1), which makes it the
#' BLUE of sigma from the two extremes (and the S1 BLUE for a median symmetric
#' about the center, whose weight is zero); under "wan_luo_shi" the Blom divisor
#' xi(n) = 2 qnorm((n - 0.375) / (n + 0.25)) of Wan et al. (2014), the S1 SD
#' estimator of that version, which also puts no weight on the median. The
#' factor C_sigma(n) = SD(X_(n) - X_(1)) / divisor, in units of sigma, comes from
#' the exact moments of the two extremes (the session's cache), so neither
#' number depends on the quantile definition `type`, which concerns a reported
#' quantile and this row has none. Returns sigma_hat and c_sigma; the caller
#' checks the inputs (finite, min < max, n >= 3).
#' @keywords internal
#' @noRd
.range_sd <- function(a, b, n, version) {
  n <- as.integer(n)
  mo <- .os_moments_cached(n, c(1L, n))
  divisor <- if (version == "blue") mo$alpha[2] - mo$alpha[1] else 2 * qnorm((n - 0.375) / (n + 0.25))
  sd_range <- sqrt(mo$B[1, 1] + mo$B[2, 2] - 2 * mo$B[1, 2])   # the SD of the range, sigma units
  list(sigma_hat = (b - a) / divisor, c_sigma = sd_range / divisor)
}
