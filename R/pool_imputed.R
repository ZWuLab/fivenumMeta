# ==============================================================================
# File: pool_imputed.R
# Purpose: Pool the studies of se_imputed_diff() in one call, under the naive and
#          under the exact standard errors side by side, by the common-effect model
#          or the DerSimonian-Laird random-effects model, with a compact print.
#
# EXPORTED FUNCTIONS:
#   - pool_imputed()        the two fits (naive, exact): estimate, SE, CI, p-value, Q, I^2, tau^2, weights
#   - print.pool_imputed()  the two-row table with its heading
#
# INTERNAL HELPERS:
#   - .pool_one()           one fit from the effects y and their variances v

#################### EXPORTED MAIN FUNCTIONS ####################

#' Pool the studies under the naive and the exact standard errors, side by side
#'
#' Takes the output of \code{\link{se_imputed_diff}} (or any data frame with the
#' columns \code{yi}, \code{sei_naive} and \code{sei}) and pools the studies
#' twice, once with the naive standard errors and once with the exact ones, by
#' the common-effect model (\code{method = "EE"}, the default) or the
#' DerSimonian-Laird random-effects model (\code{method = "DL"}). The two fits
#' share the effect estimates and differ in the variances alone, so their
#' difference shows what the exact standard errors change: the pooled estimate,
#' its standard error and confidence interval, the heterogeneity statistic
#' \eqn{\mathcal Q} and the index \eqn{I^2}, and the weights of the studies.
#'
#' @param data A data frame with the columns \code{yi} (the effect of each
#'   study), \code{sei_naive} (its standard error from the naive standard errors)
#'   and \code{sei} (from the exact ones), as \code{\link{se_imputed_diff}}
#'   returns it; the three must be numeric. A row whose effect is not finite, or
#'   whose standard error in either column is not finite and positive, is dropped
#'   from both fits, so that the two rows compare the same studies, with a
#'   warning that names the row and the column.
#' @param method \code{"EE"} (default), the common-effect (equal-effects) model,
#'   or \code{"DL"}, the random-effects model with the DerSimonian-Laird
#'   estimator of \eqn{\tau^2}. The codes are those of \code{metafor::rma()}.
#' @param level The confidence level of the interval, a number strictly between
#'   0 and 1; 0.95 by default.
#' @param x An object returned by \code{pool_imputed()}.
#' @param digits The number of significant digits the print shows.
#' @param ... Passed to \code{\link{print.data.frame}} for the table.
#' @return An object of class \code{"pool_imputed"}: a list, used as any list,
#'   whose print shows the heading and the two-row table. Its components:
#'   \code{fits}, a data frame with the rows \code{naive} and \code{exact} and
#'   the columns \code{estimate}, \code{se}, \code{ci.lb}, \code{ci.ub},
#'   \code{pval}, \code{Q}, \code{I2} (in percent) and \code{tau2}; \code{weights},
#'   a data frame with one row per study used and the columns \code{naive} and
#'   \code{exact}, the weights in percent of their total; \code{k}, the number of
#'   studies used; \code{method}; \code{level}; and \code{measure}, read from the
#'   attribute that \code{\link{se_imputed_diff}} sets (\code{"MD"} or
#'   \code{"lnVR"}), \code{NA} when \code{data} carries none (a column subset
#'   \code{data[, cols]} drops the attribute, a row subset keeps it; without it
#'   the print says "Pooled estimate").
#' @details
#' With \eqn{v_i} the squared standard error of study \eqn{i}, \eqn{w_i = 1/v_i}
#' and \eqn{K} studies, the common-effect model gives
#' \deqn{\hat\theta = \frac{\sum_i w_i y_i}{\sum_i w_i}, \qquad
#'   \mathrm{SE}(\hat\theta) = \Big(\sum_i w_i\Big)^{-1/2},}
#' \deqn{\mathcal Q = \sum_i w_i (y_i - \hat\theta)^2, \qquad
#'   I^2 = \max\Big\{0, \frac{\mathcal Q - (K - 1)}{\mathcal Q}\Big\},}
#' with \eqn{\tau^2 = 0}. \eqn{\mathcal Q} and \eqn{I^2} are the heterogeneity
#' measures of Higgins and Thompson (2002): when the studies share one effect and
#' the \eqn{v_i} are the variances of the \eqn{y_i}, \eqn{\mathcal Q} has
#' expectation \eqn{K - 1}, so \eqn{I^2} is the share of \eqn{\mathcal Q} beyond
#' that. The DerSimonian-Laird model (DerSimonian and Laird, 1986) keeps
#' \eqn{\mathcal Q} and \eqn{I^2} from these weights and sets
#' \deqn{\hat\tau^2 = \max\Big\{0, \frac{\mathcal Q - (K - 1)}{\sum_i w_i - \sum_i w_i^2 / \sum_i w_i}\Big\},}
#' \deqn{w_i^\ast = \frac{1}{v_i + \hat\tau^2},}
#' and takes \eqn{\hat\theta} and \eqn{\mathrm{SE}(\hat\theta)} as above with
#' \eqn{w_i^\ast} in place of \eqn{w_i}. The confidence interval is
#' \eqn{\hat\theta \pm z_{1 - \alpha/2}\,\mathrm{SE}(\hat\theta)} and the
#' p-value is the two-sided normal test of \eqn{\theta = 0}, as the default of
#' \code{metafor::rma()}. The weights are \eqn{w_i} under \code{"EE"} and
#' \eqn{w_i^\ast} under \code{"DL"}, in percent of their total.
#'
#' Everything is computed twice, with \eqn{v_i} the square of \code{sei_naive}
#' and with \eqn{v_i} the square of \code{sei}. Under the naive standard errors
#' an imputed mean or SD is treated as a sample mean or SD, so the variances of
#' the studies that report a summary are too small; the exact row shows the fit
#' with the variances of Wu, Shi and Tong (2026). Under the common-effect model
#' the pooled standard error is \eqn{(\sum_i w_i)^{-1/2}} and \eqn{\mathcal Q}
#' is the minimum over \eqn{\theta} of \eqn{\sum_i w_i (y_i - \theta)^2}, so
#' when no exact variance is below its naive one (the case for the difference of
#' means, since \eqn{\kappa_\mu(n) \ge 1}) the exact row has the larger pooled
#' standard error and the smaller \eqn{\mathcal Q} and \eqn{I^2}. Under
#' \code{"DL"} the direction of the pooled standard error is not fixed, because
#' the re-estimated \eqn{\hat\tau^2} enters every weight and a smaller
#' \eqn{\hat\tau^2} can outweigh the larger \eqn{v_i}.
#'
#' What this function does not supply: a p-value for \eqn{\mathcal Q}, REML and
#' the other estimators of \eqn{\tau^2}, the Knapp-Hartung adjustment and
#' moderators, which \code{metafor::rma(yi, sei = sei, data = d)} fits from the
#' same columns (Viechtbauer, 2010). The exact reference distribution of the
#' pooled estimate under imputed standard errors is the subject of a later
#' version: the intervals and p-values here are normal approximations, as
#' metafor's default ones are.
#' @examples
#' # a meta-analysis of means: controls minus cases, the data set vitd
#' d <- se_imputed_diff(vitd, suffixes = c(".controls", ".cases"))
#' fit <- pool_imputed(d)
#' fit
#' fit$fits["exact", c("estimate", "se")]
#' round(fit$weights, 1)                      # each study's weight, percent
#' pool_imputed(d, method = "DL", level = 0.9)
#' \dontrun{
#' # a meta-analysis of SDs: the log variability ratio on metamedian::dat.age
#' a <- se_imputed_diff(metamedian::dat.age, measure = "lnVR")
#' pool_imputed(a)
#' metafor::rma(yi, sei = sei, data = a, method = "EE")   # the same exact fit
#' }
#' @seealso \code{\link{se_imputed_diff}}, which builds \code{data};
#'   \code{\link{se_imputed_df}}.
#' @references
#' DerSimonian R, Laird N (1986). Meta-analysis in clinical trials. Controlled
#' Clinical Trials 7(3), 177-188.
#'
#' Higgins JPT, Thompson SG (2002). Quantifying heterogeneity in a meta-analysis.
#' Statistics in Medicine 21(11), 1539-1558.
#'
#' Viechtbauer W (2010). Conducting meta-analyses in R with the metafor package.
#' Journal of Statistical Software 36(3), 1-48.
#'
#' Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard
#' deviation from a five-number summary. Manuscript.
#' @importFrom stats pnorm qnorm
#' @export
pool_imputed <- function(data, method = c("EE", "DL"), level = 0.95) {
  method <- match.arg(method)
  if (!is.data.frame(data))
    stop("data must be a data frame with the columns yi, sei_naive and sei, as se_imputed_diff() returns it")
  need <- c("yi", "sei_naive", "sei")
  miss <- setdiff(need, names(data))
  if (length(miss))
    stop("data lacks the column", if (length(miss) > 1) "s " else " ", paste(miss, collapse = ", "),
         "; se_imputed_diff() returns them")
  for (nm in need)   # a factor or a character column would be read as codes or coerced with NAs
    if (!is.numeric(data[[nm]]))
      stop("column ", nm, " must be numeric, not ", class(data[[nm]])[1])
  if (!is.numeric(level) || length(level) != 1 || !is.finite(level) || level <= 0 || level >= 1)
    stop("level must be one number strictly between 0 and 1, such as 0.95")
  y <- as.numeric(data$yi); vn <- as.numeric(data$sei_naive)^2; ve <- as.numeric(data$sei)^2
  # a usable row has a finite effect and standard errors whose squares are positive and finite
  # (a zero would be an infinite weight); the first failing column names the dropped row
  fail <- ifelse(!is.finite(y), "yi",
          ifelse(!is.finite(vn) | vn <= 0, "sei_naive",
          ifelse(!is.finite(ve) | ve <= 0, "sei", NA_character_)))
  ok <- is.na(fail)
  if (any(!ok))
    warning("rows dropped, ", sum(!ok), " of ", length(ok), ", for a non-finite effect or a ",
            "standard error that is not finite and positive: ",
            paste0(which(!ok), " (", fail[!ok], ")", collapse = ", "))
  k <- sum(ok)
  if (k < 2) stop("pool_imputed() needs at least two usable studies; got ", k)
  y <- y[ok]; vn <- vn[ok]; ve <- ve[ok]
  fn <- .pool_one(y, vn, method, level)   # the naive variances
  fe <- .pool_one(y, ve, method, level)   # the exact variances
  fits <- as.data.frame(rbind(naive = fn$fit, exact = fe$fit))
  weights <- data.frame(naive = fn$w, exact = fe$w)
  rn <- rownames(data)[ok]
  if (!is.null(rn)) rownames(weights) <- rn
  measure <- attr(data, "measure")
  if (is.null(measure)) measure <- NA_character_
  structure(list(fits = fits, weights = weights, k = k, method = method, level = level,
                 measure = measure),
            class = c("pool_imputed", "list"))
}

#' @rdname pool_imputed
#' @export
print.pool_imputed <- function(x, digits = 4L, ...) {
  what <- switch(as.character(x$measure),
                 MD = "Pooled difference of means (group 1 minus group 2)",
                 lnVR = "Pooled log variability ratio (SD of group 1 over SD of group 2)",
                 "Pooled estimate")
  model <- switch(x$method, EE = "common-effect model (EE)",
                  DL = "random-effects model, DerSimonian-Laird (DL)")
  cat(what, ", ", model, ", ", x$k, " studies\n", sep = "")
  tab <- x$fits
  if (x$method == "EE") tab$tau2 <- NULL   # tau^2 is zero by the model, not estimated
  rownames(tab) <- c("naive SE", "exact SE")
  print(tab, digits = digits, ...)
  cat(format(100 * x$level), "% confidence intervals and p-values from the normal distribution. ",
      "Q has ", x$k - 1, " degrees of freedom, and I2 is in percent.\n", sep = "")
  if (identical(x$measure, "lnVR"))
    cat("exp() of the estimate and of the limits gives the ratio of SDs.\n")
  cat("The weights in percent are $weights. Every field with str() or unclass().\n")
  invisible(x)
}

#################### INTERNAL HELPER FUNCTIONS ####################

#' One pooled fit from the effects y and their variances v
#'
#' Returns a list with \code{fit}, the named vector (estimate, se, ci.lb, ci.ub,
#' pval, Q, I2, tau2), and \code{w}, the weights in percent. Under "EE" the
#' weights are 1/v; under "DL" they are 1/(v + tau2) with the DerSimonian-Laird
#' tau2, while Q and I2 come from the weights 1/v in both cases.
#' @keywords internal
#' @noRd
.pool_one <- function(y, v, method, level) {
  # the weights relative to the largest, c0 / v in (0, 1], with c0 the smallest variance: no
  # square overflows or underflows whatever the unit of y, and c0 restores the scale below
  c0 <- min(v)
  wt <- c0 / v; sw <- sum(wt)
  est_fe <- sum(wt * y) / sw                 # the common-effect estimate (the scale of wt cancels)
  k <- length(y); df <- k - 1
  Q <- sum(wt * (y - est_fe)^2) / c0         # Cochran's Q with the weights 1/v
  I2 <- if (Q > df) 100 * (Q - df) / Q else 0   # Higgins and Thompson, in percent (0 when Q <= K - 1)
  if (method == "EE") {
    tau2 <- 0; wst <- wt
  } else {
    # DerSimonian and Laird: tau^2 = (Q - df) / D with D = sum w - sum w^2 / sum w, truncated at
    # zero. With u = w / sum w, D = (sum w) (1 - sum u^2) = 2 (sum w) sum_{i<j} u_i u_j, and the
    # last form adds positive terms where the first cancels when one study carries the weight
    u <- wt / sw
    cross <- sum(u * (cumsum(u) - u))        # sum over i < j of u_i u_j
    D <- 2 * (sw / c0) * cross
    tau2 <- max(0, (Q - df) / D)
    wst <- c0 / (v + tau2)                   # the random-effects weights, relative as above
  }
  est <- sum(wst * y) / sum(wst)
  se <- sqrt(c0 / sum(wst))                  # 1 / sqrt(sum of 1 / (v + tau^2))
  if (!is.finite(Q) || !is.finite(se) || !is.finite(est))
    stop("the pooled fit is not representable in double precision (the effects or the variances ",
         "are too large, too small or too far apart); rescale yi and the standard errors")
  z <- qnorm((1 - level) / 2, lower.tail = FALSE)   # the upper tail directly, finite near level = 1
  fit <- c(estimate = est, se = se, ci.lb = est - z * se, ci.ub = est + z * se,
           pval = 2 * pnorm(-abs(est / se)), Q = Q, I2 = I2, tau2 = tau2)
  list(fit = fit, w = 100 * wst / sum(wst))
}
