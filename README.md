
# fivenumMeta

Many studies in a meta-analysis report a median with the quartiles or the extremes instead of
a mean and a standard deviation (SD). The established estimators (Wan, Luo and Shi, and the best
linear unbiased estimator, the BLUE, from the reported values) impute a mean and an SD from such
a summary, and the pooling then treats them as a sample mean and SD, with the standard error
sd / sqrt(n). An imputed mean is a combination of three or five reported values and is less
precise than a sample mean. Under normality, for the BLUE with the quartiles by R's default
definition and sample sizes from 5 to 1001, its standard deviation is 6% to 8% larger
than sigma / sqrt(n) when the three quartiles are reported, and that of the imputed SD is 65% to
68% larger than a sample SD's. When the minimum, median and maximum are reported instead,
both excesses grow with the sample size: from 3% to 24% for the mean, and from 9%
to more than threefold for the SD. With all five values reported the excesses are smaller, at
most 7% for the mean and 6% to 51% for the SD.

fivenumMeta computes those standard deviations exactly, from the means and covariances of the
normal order statistics that the package `osmoments` computes, returns the plug-in standard
errors that belong in the weights of a meta-analysis, and pools the studies under the naive and
the exact standard errors side by side. It covers the three data-reporting
scenarios (the minimum, median and maximum, the three quartiles, or all five), the BLUE (the default)
and the estimators of Wan, Luo and Shi, at every sample size, and it gives closed-form shortcut
formulas for the factors, for a reader without the package and for comparison.

## Installation

```r
# install.packages("remotes")
remotes::install_github("ZWuLab/osmoments")      # the engine, once it is public
remotes::install_github("ZWuLab/fivenumMeta")    # this package
```

## Usage

A meta-analysis data set with two groups per row, in the column layout of the data sets that
the `metamedian` package ships (`n.g1`, `q1.g1`, `med.g1`, `q3.g1`, `mean.g1`, `sd.g1`, and the
same with `.g2`). The example needs the packages `metamedian` (the data set) and `metafor` (the
pooling), which fivenumMeta itself does not require:

```r
library(fivenumMeta)
d <- se_imputed_diff(
  data     = metamedian::dat.age,   # one row per study; a study reports quartiles or a mean and SD
  suffixes = c(".g1", ".g2"),       # the two groups' column suffixes, group 1 first
  version  = "blue"                 # the estimator version: the BLUE (default) or "wan_luo_shi"
)
pool_imputed(d)                                       # the pooled difference under the naive and the exact SEs
metafor::rma(yi, sei = sei, data = d, method = "EE")  # the same exact fit with metafor
metafor::rma(yi, sei = sei, data = d)                 # random effects (REML) with metafor, which pool_imputed() does not fit
```

`pool_imputed()` fits the common-effect model by default and the DerSimonian-Laird
random-effects model with `method = "DL"`. For a meta-analysis of SDs,
`se_imputed_diff(..., measure = "lnVR")` gives the log variability ratio of the two groups
with its two standard errors. For one group, `se_imputed_df()` appends to each row the imputed
mean and SD, the naive standard error, the exact plug-in standard errors of the mean and of
the SD, and the scenario it read off the row, a study's sample mean and SD kept as they are
and a mean with the range given the SD imputed from the range. For one study,
`se_imputed(c(q1 = 10, med = 12, q3 = 15), n = 25, scenario = "S2")`. The factors themselves
are `c_mu()`, `c_sigma()` and `inflation_factors()`, and the shortcut formulas
`c_mu_shortcut()` and `c_sigma_shortcut()` need no integration.

The factors are exact at every sample size. At n = 4Q + 1 the reported quartiles are single
order statistics. At any other n each is a fixed combination of two neighboring order
statistics, by R's default definition (type 7 of `quantile()`) unless the argument `type`
names another, and the factors there depend on that definition.

The vignette, `vignette("fivenumMeta")`, walks through the example above and explains the
columns, the scenarios and the convention on the sample size. A second vignette,
`vignette("two-meta-analyses")`, reproduces the two meta-analyses of Wu, Shi and Tong (2026), one
of means on the data set `vitd` that ships with the package and one of SDs on
`metamedian::dat.age`, each in three lines.

## The papers

Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard deviation from a
five-number summary. Manuscript.

It defines the factors, gives their limits in the sample size, fits the shortcut formulas and
reports the checks of the computation. `citation("fivenumMeta")` prints this reference. The
exact moments of the order statistics come from the package `osmoments`, described in:

Wu Z, Tong T (2026). osmoments: Exact means, variances and covariances of order statistics in
R. Manuscript.

> Developed with AI assistance; see [AI-USE.md](AI-USE.md).
