# fivenumMeta: Five-Number Summaries in Meta-Analysis, Imputed Means and SDs with Standard Errors

Meta-analyses often pool studies that report a median with quartiles or
extremes rather than a mean and a standard deviation, and the
established estimators impute the mean and the SD from those summaries.
An imputed mean is not a sample mean, and its standard deviation is not
\\\sigma/\sqrt{n}\\. This package computes, under normality, the
standard deviation factors \\C\_\mu(n)\\ and \\C\_\sigma(n)\\ such that
\\SD(\hat\mu) = C\_\mu(n)\\\sigma\\ and \\SD(\hat\sigma) =
C\_\sigma(n)\\\sigma\\, and the inflation factors \\\kappa\_\mu(n) =
C\_\mu(n)\sqrt{n}\\ and \\\kappa\_\sigma(n) = C\_\sigma(n)/\sqrt{1 -
b_S(n)^2}\\, exactly, from the means and covariances of the normal order
statistics that the package osmoments computes, at every sample size.
The estimators covered are the best linear unbiased estimator (BLUE) of
the reported values, the default, and the estimators of Wan, Luo and
Shi, under the data-reporting scenarios `"S1"` (minimum, median,
maximum), `"S2"` (the three quartiles), `"S3"` (all five), and, for the
BLUE, `"S4"` and `"S5"`, which add the 10th and 90th percentiles to
`"S2"` and `"S3"`.

## Details

For a meta-analysis data set, one row per study,
[`se_imputed_df`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_df.md)
reads the scenario off the columns each row has and appends the imputed
mean and SD with their naive and exact plug-in standard errors, and
[`se_imputed_diff`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_diff.md)
does the same for two groups per row and adds the effect of the study,
the difference of means or the log variability ratio, with its naive and
its exact standard error, in the names that
[`metafor::rma()`](https://wviechtb.github.io/metafor/reference/rma.uni.html)
takes as they are.
[`pool_imputed`](https://zwulab.github.io/fivenumMeta/reference/pool_imputed.md)
pools those studies in one call, under the naive and the exact standard
errors side by side, by the common-effect or the DerSimonian-Laird
model. For one study,
[`se_imputed`](https://zwulab.github.io/fivenumMeta/reference/se_imputed.md)
returns the same quantities from the reported values and the sample
size. The factors themselves are
[`c_mu`](https://zwulab.github.io/fivenumMeta/reference/c_mu.md),
[`c_sigma`](https://zwulab.github.io/fivenumMeta/reference/c_sigma.md)
and
[`inflation_factors`](https://zwulab.github.io/fivenumMeta/reference/inflation_factors.md),
with the closed-form shortcut formulas
[`c_mu_shortcut`](https://zwulab.github.io/fivenumMeta/reference/c_mu_shortcut.md)
and
[`c_sigma_shortcut`](https://zwulab.github.io/fivenumMeta/reference/c_sigma_shortcut.md)
for use without the integration. The weights and the reported ranks are
[`estimator_weights`](https://zwulab.github.io/fivenumMeta/reference/estimator_weights.md)
and
[`scenario_ranks`](https://zwulab.github.io/fivenumMeta/reference/scenario_ranks.md),
and
[`imputed_moments`](https://zwulab.github.io/fivenumMeta/reference/imputed_moments.md)
gives the variances and the covariance of the two imputed statistics.
Every function defaults to the BLUE; the estimators of Wan, Luo and Shi
are reached by `version = "wan_luo_shi"` in the study-level functions
and in
[`inflation_factors()`](https://zwulab.github.io/fivenumMeta/reference/inflation_factors.md),
and by the argument `method` in the others.

## References

Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard
deviation from a five-number summary. Manuscript. The paper this package
accompanies: the factors, their limits, the shortcut formulas and the
checks of the computation.

Wu Z (2026). osmoments: Exact Moments of Order Statistics. R package
version 0.1.0. The engine that computes the moments of the normal order
statistics.

Wu Z, Tong T (2026). osmoments: Exact means, variances and covariances
of order statistics in R. Manuscript. The paper of the engine.

Wan X, Wang W, Liu J, Tong T (2014). Estimating the sample mean and
standard deviation from the sample size, median, range and/or
interquartile range. BMC Medical Research Methodology 14, 135.

Luo D, Wan X, Liu J, Tong T (2018). Optimally estimating the sample mean
from the sample size, median, mid-range, and/or mid-quartile range.
Statistical Methods in Medical Research 27(6), 1785-1805.

Shi J, Luo D, Weng H, Zeng XT, Lin L, Chu H, Tong T (2020). Optimally
estimating the sample standard deviation from the five-number summary.
Research Synthesis Methods 11(5), 641-654.

Balakrishnan N, Rychtar J, Taylor D, Walter SD (2022). Unified approach
to optimal estimation of mean and standard deviation from sample
summaries. Statistical Methods in Medical Research 31(11), 2087-2103.

McGrath S, Katzenschlager S, Zimmer AJ, Seitel A, Steele R, Benedetti A
(2023). Standard error estimation in meta-analysis of studies reporting
medians. Statistical Methods in Medical Research 32(2), 373-388.

## Author

**Maintainer**: Zheyang Wu <zheyangwu@wpi.edu>
