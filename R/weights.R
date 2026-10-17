# ==============================================================================
# File: weights.R
# Purpose: The weight vectors of the named estimators (Wan, Luo, Shi and the BLUE) and
#          the order-statistic ranks reported under the scenarios S1 to S5.
#
# EXPORTED FUNCTIONS:
#   - estimator_weights() weight vectors of the named estimators (wan/luo/shi/blue)
#   - scenario_ranks()    order-statistic ranks reported under S1 to S5 at the exact-rank sizes
#
# INTERNAL HELPERS:
#   - .reported_design()   the reported values as linear combinations of order statistics
#   - .reported_moments()  their means and covariance matrix (standard normal units)
#   - .check_moments()     a user-supplied `moments` object against the design
#   - .os_moments_cached() the engine's moments at (n, ranks), kept for the session
#   - .blue_weights()      Lloyd (1952) GLS weights from (alpha, B)

#################### EXPORTED MAIN FUNCTIONS ####################

#' Weight vectors of the published quantile-based estimators and the BLUE
#'
#' Returns the weight vector \eqn{w} such that the estimator is \eqn{w'Y}, with
#' \eqn{Y} the reported values in increasing order (for example the minimum, the
#' median and the maximum under S1). All formulas were read from the source PDFs.
#'
#' @inheritParams imputed_moments
#' @param n A single integer sample size. The Wan, Luo and Shi presets are closed
#'   forms in \eqn{n} and accept any \eqn{n \ge 2}; the BLUE needs the sample
#'   sizes of \code{\link{imputed_moments}}.
#' @param estimand \code{"mean"} or \code{"sd"}.
#' @param method One of the following; \code{"blue"} is the default. Under S4
#'   and S5 only \code{"blue"}.
#'   \describe{
#'   \item{\code{"blue"}}{The default. Lloyd's generalized-least-squares best
#'     linear unbiased estimator computed from the exact means and covariance
#'     matrix of the reported values (Balakrishnan 2022 eqs. 8-9; Yang 2022
#'     Definition 3.1).}
#'   \item{\code{"wan"}}{The estimators used/recommended in Wan et al. (2014).
#'     Mean: Hozo's \eqn{(a+2m+b)/4} for S1 (Wan eq. 3), \eqn{(q_1+m+q_3)/3} for
#'     S2 (Wan eq. 14), Bland's \eqn{(a+2q_1+2m+2q_3+b)/8} for S3 (Wan eq. 10).
#'     SD: \eqn{(b-a)/\xi_B(n)} for S1 (eq. 9), \eqn{(q_3-q_1)/\eta_B(n)} for S2
#'     (eq. 16), the average of the two for S3 (eq. 13), with Blom-approximated
#'     divisors \eqn{\xi_B(n) = 2\Phi^{-1}\{(n-0.375)/(n+0.25)\}},
#'     \eqn{\eta_B(n) = 2\Phi^{-1}\{(0.75n-0.125)/(n+0.25)\}}.}
#'   \item{\code{"luo"}}{Mean only. Luo et al. (2018) eqs. (7), (11), (15):
#'     S1 \eqn{w_1 = 4/(4+n^{0.75})} on the mid-range; S2 \eqn{w = 0.7+0.39/n} on
#'     the mid-quartile range; S3 \eqn{w_1 = 2.2/(2.2+n^{0.75})},
#'     \eqn{w_2 = 0.7-0.72/n^{0.55}}.}
#'   \item{\code{"shi"}}{SD only, S3 only. Shi et al. (2020) eqs. (5), (9), (10):
#'     weight \eqn{w = 1/(1+0.07 n^{0.6})} on \eqn{(b-a)/\xi_B}, rest on
#'     \eqn{(q_3-q_1)/\eta_B}. For S1/S2 Shi recommends Wan's estimators.}
#'   }
#'
#' @return Named numeric weight vector over the reported values, in increasing
#'   order, named as the scenario reports them: \code{a} the minimum, \code{p10}
#'   the 10th percentile, \code{q1}, \code{m} the median, \code{q3}, \code{p90}
#'   the 90th percentile and \code{b} the maximum.
#' @details The Wan/Luo/Shi weights are closed-form functions of \eqn{n} and
#'   are returned for any single integer \eqn{n \ge 2}, as the source papers
#'   apply them. The \code{"blue"} weights are computed from the exact means and
#'   covariance matrix of the reported values at any \eqn{n} (see
#'   \code{\link{imputed_moments}}, argument \code{type}) and use \code{moments}.
#' @examples
#' estimator_weights(25, "S3", "mean")           # the BLUE of the mean, the default
#' estimator_weights(25, "S3", "sd")             # the BLUE of the SD
#' estimator_weights(25, "S3", "mean", "luo")
#' estimator_weights(25, "S3", "sd", "shi")
#' @references
#' Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard
#' deviation from a five-number summary. Manuscript.
#' @export
estimator_weights <- function(n, scenario = c("S3", "S1", "S2", "S4", "S5"),
                        estimand = c("mean", "sd"),
                        method = c("blue", "luo", "wan", "shi"),
                        moments = NULL, type = 7) {
  scenario <- match.arg(scenario); estimand <- match.arg(estimand)
  method <- match.arg(method)
  # the BLUE needs the moments of the reported values; the formula presets only need an n
  if (method == "blue") .check_scenario_n(n, scenario) else {
    .check_n_integer(n)
    if (n < 2) stop("n must be a single integer >= 2; got: ", n)
  }

  nm <- switch(scenario, S1 = c("a", "m", "b"), S2 = c("q1", "m", "q3"),
               S3 = c("a", "q1", "m", "q3", "b"),
               S4 = c("p10", "q1", "m", "q3", "p90"),
               S5 = c("a", "p10", "q1", "m", "q3", "p90", "b"))

  if (method == "blue") {
    rm <- .reported_moments(n, scenario, type, moments)
    if (!rm$full_rank)
      stop("the reported values are linearly dependent at n = ", n, " under ", scenario,
           " and quantile type ", type, " (two of them are the same order statistic), ",
           "so the BLUE is not defined")
    bw <- .blue_weights(rm$alpha, rm$B)
    w <- if (estimand == "mean") bw$w_mu else bw$w_sig
    names(w) <- nm
    return(w)
  }
  if (scenario %in% c("S4", "S5"))
    stop("the Wan-Luo-Shi estimators are not defined for ", scenario,
         "; use method = \"blue\"")

  # Divisors shared by the wan/shi SD presets: Blom's (1958) plotting-position
  # approximation E[Z_(r:n)] ~ qnorm((r - 0.375)/(n + 0.25)) applied to the
  # extreme and upper-quartile ranks (Wan 2014 eq. 8; bib key blom1958statistical)
  xi_b  <- 2 * qnorm((n - 0.375) / (n + 0.25))         # ~ 2 E[Z_(n:n)],   Wan eq. 9
  eta_b <- 2 * qnorm((0.75 * n - 0.125) / (n + 0.25))  # ~ 2 E[Z_(3Q+1:n)], Wan eq. 16

  w <- switch(paste(estimand, method, sep = "."),
    mean.wan = switch(scenario,
      S1 = c(1, 2, 1) / 4,                # Hozo (a+2m+b)/4, Wan eq. 3
      S2 = c(1, 1, 1) / 3,                # Wan eq. 14
      S3 = c(1, 2, 2, 2, 1) / 8),         # Bland, Wan eq. 10
    mean.luo = switch(scenario,
      S1 = { w1 <- 4 / (4 + n^0.75);  c(w1 / 2, 1 - w1, w1 / 2) },       # Luo eq. 7
      S2 = { w1 <- 0.7 + 0.39 / n;    c(w1 / 2, 1 - w1, w1 / 2) },       # Luo eq. 11
      S3 = { w1 <- 2.2 / (2.2 + n^0.75); w2 <- 0.7 - 0.72 / n^0.55      # Luo eq. 15
             c(w1 / 2, w2 / 2, 1 - w1 - w2, w2 / 2, w1 / 2) }),
    sd.wan = switch(scenario,
      S1 = c(-1, 0, 1) / xi_b,                                            # Wan eq. 9
      S2 = c(-1, 0, 1) / eta_b,                                           # Wan eq. 16
      S3 = c(-1 / (2 * xi_b), -1 / (2 * eta_b), 0,                        # Wan eq. 13
              1 / (2 * eta_b),  1 / (2 * xi_b))),
    sd.shi = switch(scenario,
      S3 = { ws <- 1 / (1 + 0.07 * n^0.6)                                 # Shi eq. 9
             c(-ws / xi_b, -(1 - ws) / eta_b, 0, (1 - ws) / eta_b, ws / xi_b) },
      stop("Shi et al. (2020) define the SD estimator for S3 only; ",
           "use method = 'wan' for S1/S2.")),
    mean.shi = stop("Shi et al. (2020) is an SD estimator; no mean preset."),
    sd.luo   = stop("Luo et al. (2018) is a mean estimator; no SD preset."))
  names(w) <- nm
  w
}

#' Ranks of the order statistics reported under each scenario
#'
#' For \eqn{n = 4Q+1}: S1 reports ranks \eqn{(1, 2Q+1, 4Q+1)} (min, median, max),
#' S2 ranks \eqn{(Q+1, 2Q+1, 3Q+1)} (the quartiles), S3 all five
#' (Luo 2018 sec. 2.2; McGrath 2023 scenario labels). For \eqn{n = 20m+1}: S4
#' reports the quartiles with the 10th and 90th percentiles at ranks
#' \eqn{2m+1} and \eqn{18m+1}, that is \eqn{p(n-1)+1} for \eqn{p = 0.1, 0.9}, and
#' S5 adds the minimum and maximum (Balakrishnan 2022 sec. 3.1, Scenarios 5 and 7).
#'
#' These are the sample sizes at which every reported value is a single order
#' statistic under R's default quantile definition (type 7). At other sample
#' sizes the functions of the package build the reported values from two
#' adjacent order statistics each (see \code{\link{imputed_moments}}).
#'
#' @inheritParams imputed_moments
#' @param n A single integer sample size at which every reported value is one
#'   order statistic under type 7: odd \eqn{n \ge 3} under S1, \eqn{n = 4Q + 1}
#'   with \eqn{Q \ge 1} under S2 and S3, and \eqn{n = 20m + 1} with \eqn{m \ge 1}
#'   under S4 and S5. Any other \eqn{n} is an error here, while the other
#'   functions of the package accept it.
#' @return Integer vector of ranks, increasing.
#' @examples scenario_ranks(25, "S3")
#' @references
#' Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard
#' deviation from a five-number summary. Manuscript.
#' @export
scenario_ranks <- function(n, scenario = c("S3", "S1", "S2", "S4", "S5")) {
  scenario <- match.arg(scenario)
  .check_exact_rank_n(n, scenario)
  if (scenario == "S1") return(as.integer(c(1, (n + 1) / 2, n)))
  Q <- (n - 1) / 4
  if (scenario == "S2") return(as.integer(c(Q + 1, 2 * Q + 1, 3 * Q + 1)))
  if (scenario == "S3") return(as.integer(c(1, Q + 1, 2 * Q + 1, 3 * Q + 1, 4 * Q + 1)))
  m <- (n - 1) / 20                      # n = 20m+1: the 10th and 90th percentiles at 2m+1, 18m+1
  core <- c(2 * m + 1, Q + 1, 2 * Q + 1, 3 * Q + 1, 18 * m + 1)
  as.integer(if (scenario == "S4") core else c(1, core, n))
}

#################### INTERNAL HELPER FUNCTIONS ####################

#' The reported values as linear combinations of order statistics
#'
#' Under scenario `scenario` at sample size `n`, each reported value is an order
#' statistic (the minimum, the maximum) or a sample quantile under definition
#' `type`, which is a fixed linear combination (1 - h) X_(j) + h X_(j+1) of two
#' adjacent order statistics (osmoments::os_quantile_weights()). Returns the
#' ranks that carry a weight and the matrix W with one column per reported value,
#' so that the reported vector is Y = W' X_(ranks). A weight below 1e-12 in
#' absolute value is a rounding residue of the quantile arithmetic and is set to
#' zero, so that its rank is not computed, and the column is renormalized to sum
#' to one.
#' @keywords internal
#' @noRd
.reported_design <- function(n, scenario, type = 7) {
  if (length(type) != 1 || !is.numeric(type) || !(type %in% 1:9))
    stop("type must be one of the quantile definitions 1 to 9 of stats::quantile(); got: ",
         paste(type, collapse = ", "))
  # the probabilities of the reported quantiles, and whether the extremes are reported
  p <- switch(scenario, S1 = 0.5, S2 = c(0.25, 0.5, 0.75), S3 = c(0.25, 0.5, 0.75),
              S4 = c(0.1, 0.25, 0.5, 0.75, 0.9), S5 = c(0.1, 0.25, 0.5, 0.75, 0.9))
  extremes <- scenario %in% c("S1", "S3", "S5")
  qw <- osmoments::os_quantile_weights(p, n, type)
  Wq <- qw$W
  Wq[abs(Wq) < 1e-12] <- 0                              # rounding residues of (j, h)
  Wq <- sweep(Wq, 2, colSums(Wq), "/")                  # each quantile's weights sum to one again
  ranks <- sort(unique(c(if (extremes) c(1, n), qw$ranks[rowSums(Wq != 0) > 0])))
  W <- matrix(0, length(ranks), length(p) + 2 * extremes)
  qcols <- if (extremes) 1 + seq_along(p) else seq_along(p)
  W[match(qw$ranks, ranks, nomatch = 0), qcols] <- Wq[qw$ranks %in% ranks, , drop = FALSE]
  if (extremes) { W[1, 1] <- 1; W[length(ranks), ncol(W)] <- 1 }
  list(ranks = as.integer(ranks), W = W)
}

#' Means and covariance matrix of the reported values, in standard normal units
#'
#' alpha_Y = W' alpha and B_Y = W' B W, from the exact moments of the normal order
#' statistics at the ranks of .reported_design() (osmoments::os_moments()), by the
#' formulas of osmoments::os_lin_moments(). `moments`, when given, must be os_moments() at those
#' ranks and at this n (.check_moments()); when NULL, the engine's moments come from the session's
#' cache (.os_moments_cached()). At n = 4Q + 1 under type 7 the ranks are those of
#' scenario_ranks(), and W is the identity, so alpha_Y and B_Y are alpha and B. The component
#' `full_rank` says whether the reported values are linearly independent, which the BLUE needs.
#' @keywords internal
#' @noRd
.reported_moments <- function(n, scenario, type = 7, moments = NULL) {
  des <- .reported_design(n, scenario, type)
  if (is.null(moments)) moments <- .os_moments_cached(n, des$ranks)
  else .check_moments(moments, n, des$ranks, scenario, type)
  # W' alpha and W' B W, the formulas of osmoments::os_lin_moments(), taken here from the
  # components alpha and B so that a list carrying only those (and ranks) is accepted too
  a <- drop(crossprod(des$W, moments$alpha))
  B <- unname(crossprod(des$W, moments$B %*% des$W))
  # two reported values that are the same order statistic (a quantile at an extreme rank
  # under some definitions at a small n) make W rank-deficient and B_Y singular: the
  # variance w'B_Y w of fixed weights is still defined, the BLUE is not, and
  # estimator_weights() refuses it
  full_rank <- qr(des$W)$rank == ncol(des$W)
  list(alpha = unname(a), B = B, ranks = des$ranks, W = des$W, full_rank = full_rank,
       moments = moments, n = n, scenario = scenario, type = type)
}

#' Check a user-supplied `moments` object against the design
#'
#' `moments` is osmoments::os_moments() for the standard normal at sample size `n` and
#' `ranks`, or a list with the components alpha and B alone (the former interface), which
#' can be checked only for its size. Every component the object carries is compared with
#' what the design needs: the ranks, the sample size, the distribution, the location and
#' scale, and the shape and finiteness of alpha and B (B symmetric).
#' @keywords internal
#' @noRd
.check_moments <- function(moments, n, ranks, scenario, type) {
  k <- length(ranks)
  where <- paste0(" (n = ", n, ", ", scenario, ", type ", type, ")")
  if (!is.list(moments) || is.null(moments$alpha) || is.null(moments$B))
    stop("'moments' must be osmoments::os_moments() at the ranks ", paste(ranks, collapse = ", "),
         where, ", or a list with its components alpha and B")
  # a list without ranks (alpha and B alone, as the former interface took) must match in size
  given <- if (is.null(moments$ranks)) { if (length(moments$alpha) == k) ranks else -1L } else
    as.integer(moments$ranks)
  if (!identical(given, ranks))
    stop("'moments' must be osmoments::os_moments() at the ranks ", paste(ranks, collapse = ", "), where)
  # the same ranks at another n are a different set of order statistics
  if (!is.null(moments$n) && moments$n != n)
    stop("'moments' was computed at n = ", moments$n, ", not at n = ", n)
  # the factors are those of the normal model
  if (!is.null(moments$dist$name) && !identical(moments$dist$name, "norm"))
    stop("'moments' must be the moments of the normal order statistics, not of dist = \"",
         moments$dist$name, "\"")
  # a moments object computed at another location or scale would shift alpha and B
  if (!is.null(moments$location) && !isTRUE(all.equal(c(moments$location, moments$scale), c(0, 1))))
    stop("'moments' must be computed at location 0 and scale 1, the standard normal")
  alpha <- moments$alpha; B <- moments$B
  if (!is.numeric(alpha) || length(alpha) != k || any(!is.finite(alpha)))
    stop("'moments$alpha' must be ", k, " finite numbers")
  if (!is.matrix(B) || !is.numeric(B) || any(dim(B) != k) || any(!is.finite(B)) ||
      !isSymmetric(unname(B), tol = 1e-10))
    stop("'moments$B' must be a finite symmetric ", k, " x ", k, " matrix")
  invisible(TRUE)
}

# The session's cache of the engine's moments, keyed on n and the ranks, so that the rows of a
# data set, or repeated calls at one sample size, integrate once. At most .moments_cache_max
# entries; an entry is a few kilobytes. A `moments` argument supplied by the caller bypasses it.
.moments_cache <- new.env(parent = emptyenv())
.moments_cache_max <- 200L

#' osmoments::os_moments(n, ranks) for the standard normal, from the session's cache
#' @keywords internal
#' @noRd
.os_moments_cached <- function(n, ranks) {
  key <- paste(n, paste(ranks, collapse = ","))
  mo <- .moments_cache[[key]]
  if (is.null(mo)) {
    mo <- osmoments::os_moments(n, ranks)
    if (length(ls(.moments_cache)) >= .moments_cache_max)
      rm(list = ls(.moments_cache)[1L], envir = .moments_cache)
    assign(key, mo, envir = .moments_cache)
  }
  mo
}

#' Lloyd (1952) GLS/BLUE weights for location and scale from (alpha, B)
#'
#' w_mu' = (aaa 1' - a1a alpha') Omega / Delta,
#' w_sig' = (a11 alpha' - a1a 1') Omega / Delta, Omega = B^{-1},
#' Delta = (1'Omega 1)(alpha'Omega alpha) - (1'Omega alpha)^2
#' (Lloyd 1952; Balakrishnan 2022 eqs. 8-10; identical to Yang 2022
#' Definition 3.1).
#' @keywords internal
#' @noRd
.blue_weights <- function(alpha, B) {
  Om <- tryCatch(solve(B), error = function(e)
    stop("the covariance matrix of the reported values is singular, so the BLUE is not defined",
         call. = FALSE))
  one <- rep(1, length(alpha))
  a11 <- as.numeric(t(one) %*% Om %*% one)
  a1a <- as.numeric(t(one) %*% Om %*% alpha)
  aaa <- as.numeric(t(alpha) %*% Om %*% alpha)
  Delta <- a11 * aaa - a1a^2
  w_mu  <- as.numeric((aaa * t(one) - a1a * t(alpha)) %*% Om) / Delta
  w_sig <- as.numeric((a11 * t(alpha) - a1a * t(one)) %*% Om) / Delta
  list(w_mu = w_mu, w_sig = w_sig)
}
