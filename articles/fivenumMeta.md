# Standard errors for means and SDs imputed from five-number summaries

## The problem, and what this package does

Many studies in a meta-analysis report a median with the quartiles, or
with the minimum and the maximum, instead of the sample mean and sample
standard deviation (SD). The established estimators (Wan, Luo and Shi,
and the best linear unbiased estimators from the reported values) impute
population mean and standard deviation parameters from such a summary.
In the pooling stage of a meta-analysis, such imputed mean and SD
estimators are often treated as the sample mean and sample SD, with the
naive standard errors of the sample mean and sample SD. However, the
imputed estimators are generally less precise than the sample
estimators. Thus their variations/uncertainties should be accounted for.

This package computes, under normality, how much less precise. The
standard deviation of the imputed mean is $`C_\mu(n)\,\sigma`$ and that
of the imputed SD is $`C_\sigma(n)\,\sigma`$, with factors that depend
on the sample size, on which values were reported, on the estimator and
on the definition of the sample quantiles, and not on the data. The
package computes the factors exactly from the means and covariances of
the normal order statistics. It returns the plug-in standard errors
$`C_\mu(n)\,\hat\sigma`$ and $`C_\sigma(n)\,\hat\sigma`$ that belong in
the weights of a meta-analysis.

At a sample size of 25, the standard deviation of the imputed mean
exceeds $`\sigma/\sqrt n`$ by 8% when the three quartiles are reported,
by 13% when the minimum, median and maximum are reported, and by 3% when
all five values are. The standard deviation of the imputed SD exceeds
that of a sample SD by 65%, 26% and 10% in the same three cases. With
only the minimum, median and maximum reported, the excess grows with the
sample size, to 24% and 243% at a sample size of 1001.

### The evidence on a public data set

The example data set comes from the meta-analysis of Katzenschlager et
al. (2021) of the age of COVID-19 patients who died against those who
survived, shipped with the `metamedian` package as `dat.age`: 51
studies, two groups each, of which 29 report the quartiles and 22 a mean
and SD. The call below imputes the mean and SD of each group where they
are missing, keeps the sample ones where they are not, and appends the
difference of means (nonsurvivors minus survivors) with its standard
error under the naive rule and under the exact factors. The quartiles
are taken under R’s default definition (type 7), the package’s default.
A study that states another definition is given its `type`.

``` r

d <- se_imputed_diff(
  data     = metamedian::dat.age,   # one row per study
  suffixes = c(".g1", ".g2"),       # the two groups' suffixes, group 1 first
  version  = "blue"                 # the BLUE (default) or "wan_luo_shi"
)
str(d)   # the input columns, each group's output columns, the difference
#> 'data.frame':    51 obs. of  32 variables:
#>  $ author         : chr  "Russo, Vincenzo (PHRS)" "Yang, Xiaobo (Lancet)" "Shi, Qiao (DC) Diabetic cohort" "Shi, Qiao (DC) Non-diabetic cohort" ...
#>  $ n.g1           : int  35 32 31 16 29 84 134 65 42 17 ...
#>  $ q1.g1          : num  NA NA 65 68 NA 64 NA 70 NA 63 ...
#>  $ med.g1         : num  NA NA 76 72 NA 71 NA 76 NA 72 ...
#>  $ q3.g1          : num  NA NA 82 81 NA 80 NA 83 NA 81 ...
#>  $ mean.g1        : num  77 64.6 NA NA 63.6 NA 68.7 NA 75.3 NA ...
#>  $ sd.g1          : num  8.3 11.2 NA NA 17.7 ...
#>  $ n.g2           : int  157 20 122 137 351 214 315 274 149 85 ...
#>  $ q1.g2          : num  NA NA 56 56 NA 37 NA 64 NA 47 ...
#>  $ med.g2         : num  NA NA 63 63 NA 51 NA 68 NA 53 ...
#>  $ q3.g2          : num  NA NA 69 70 NA 63 NA 74 NA 66 ...
#>  $ mean.g2        : num  65.6 51.9 NA NA 52.9 ...
#>  $ sd.g2          : num  15.6 12.9 NA NA 16.4 ...
#>  $ scenario.g1    : chr  "sample" "sample" "S2" "S2" ...
#>  $ mu_hat.g1      : num  77 64.6 74.1 73.9 63.6 ...
#>  $ sigma_hat.g1   : num  8.3 11.2 13.2 10.5 17.7 ...
#>  $ se_naive.g1    : num  1.4 1.98 2.37 2.63 3.28 ...
#>  $ se_mu.g1       : num  1.4 1.98 2.53 2.79 3.28 ...
#>  $ se_sigma.g1    : num  1 1.42 2.71 3 2.35 ...
#>  $ version_used.g1: chr  "sample" "sample" "blue" "blue" ...
#>  $ exact.g1       : logi  NA NA TRUE TRUE NA TRUE ...
#>  $ scenario.g2    : chr  "sample" "sample" "S2" "S2" ...
#>  $ mu_hat.g2      : num  65.6 51.9 62.6 63 52.9 ...
#>  $ sigma_hat.g2   : num  15.6 12.9 9.75 10.49 16.35 ...
#>  $ se_naive.g2    : num  1.245 2.885 0.883 0.896 0.873 ...
#>  $ se_mu.g2       : num  1.245 2.885 0.949 0.965 0.873 ...
#>  $ se_sigma.g2    : num  0.882 2.078 1.027 1.048 0.618 ...
#>  $ version_used.g2: chr  "sample" "sample" "blue" "blue" ...
#>  $ exact.g2       : logi  NA NA TRUE TRUE NA TRUE ...
#>  $ yi             : num  11.4 12.7 11.5 10.9 10.8 ...
#>  $ sei_naive      : num  1.88 3.5 2.53 2.78 3.4 ...
#>  $ sei            : num  1.88 3.5 2.7 2.95 3.4 ...
#>  - attr(*, "measure")= chr "MD"
table(d$scenario.g1)
#> 
#>     S2 sample 
#>     29     22
range(d$se_mu.g1[d$scenario.g1 == "S2"] / d$se_naive.g1[d$scenario.g1 == "S2"])
#> [1] 1.05489 1.07746
```

The result keeps every input column and appends, for each group with its
suffix, the scenario read off the row, the imputed or sample mean and
SD, and the three standard errors, then the difference `yi` with its
standard errors `sei_naive` and `sei`. The section “Where to look
afterwards” defines each column.
[`pool_imputed()`](https://zwulab.github.io/fivenumMeta/reference/pool_imputed.md)
pools the studies under the naive and the exact standard errors side by
side, with the common-effect model by default:

``` r

fit <- pool_imputed(d)
fit
#> Pooled difference of means (group 1 minus group 2), common-effect model (EE), 51 studies
#>          estimate     se ci.lb ci.ub pval     Q    I2
#> naive SE    11.89 0.2715 11.35 12.42    0 443.7 88.73
#> exact SE    11.92 0.2853 11.36 12.48    0 401.1 87.53
#> 95% confidence intervals and p-values from the normal distribution. Q has 50 degrees of freedom, and I2 is in percent.
#> The weights in percent are $weights. Every field with str() or unclass().
```

The estimates, standard errors and interval limits are in years and
$`I^2`$ is in percent. For the random-effects model, `metafor` takes the
same columns (here with REML):

``` r

re_naive <- metafor::rma(yi, sei = sei_naive, data = d)   # random effects, REML
re_exact <- metafor::rma(yi, sei = sei, data = d)
row <- function(fit) {
  c(estimate = as.numeric(fit$b), se = fit$se, ci.lb = fit$ci.lb,
    ci.ub = fit$ci.ub, tau2 = fit$tau2, I2 = fit$I2)
}
round(rbind(`naive SE` = row(re_naive), `exact SE` = row(re_exact)), 3)
#>          estimate    se  ci.lb  ci.ub   tau2     I2
#> naive SE   13.342 0.863 11.650 15.034 30.663 88.830
#> exact SE   13.343 0.862 11.653 15.033 30.087 87.626
```

$`\tau^2`$ is in years squared.

Three things can be read off the output. First, for every study that
reports quartiles, the exact standard error of its imputed mean is 5.5%
to 7.7% larger than the naive one, so such a study carries 2% to 5% less
weight. Under the common-effect model, the share of the weight on the
quartile studies falls from 69.9% to 66.7%, and the pooled standard
error rises by 5.1%. Second, under the random-effects model the pooled
difference moves by 0.1% of its standard error, and the standard error
itself by 0.1%, because the between-study variance dominates the total
there ($`I^2`$ of 89%). Third, $`\tau^2`$ is smaller under the exact
standard errors, 30.1 against 30.7, because larger within-study
variances leave less of the spread of the estimates to be read as
heterogeneity.

The calls to use:
[`se_imputed_diff()`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_diff.md)
for two groups per row,
[`se_imputed_df()`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_df.md)
for one, and
[`pool_imputed()`](https://zwulab.github.io/fivenumMeta/reference/pool_imputed.md)
to pool, with `version = "blue"` (the default).
`version = "wan_luo_shi"` gives the estimators of Wan, Luo and Shi
instead, for comparability with an earlier analysis. For a meta-analysis
of SDs, `se_imputed_diff(..., measure = "lnVR")` gives the log
variability ratio of the two groups with its two standard errors, as the
second vignette shows.

## Setup: the columns and the scenarios

[`se_imputed_df()`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_df.md)
takes one row per study and reads these columns, each named by an
argument with the default shown, with a `suffix` appended for data sets
that carry two groups side by side:

| Argument | Default column | Holds |
|----|----|----|
| `n` | `n` | the sample size |
| `min`, `max` | `min`, `max` | the minimum and the maximum |
| `q1`, `med`, `q3` | `q1`, `med`, `q3` | the first quartile, the median, the third quartile |
| `mean`, `sd` | `mean`, `sd` | the sample mean and SD, where a study gives them |

A column that the data set does not have, or an argument set to `NULL`,
counts as absent for every row, and a missing value counts as absent for
that row. The scenario of a row is read off the values present: all five
give $`S_3`$; the three quartiles give $`S_2`$; the minimum, median and
maximum give $`S_1`$; a sample mean and SD with no summary give
“sample”, whose standard error is sd / sqrt(n); and a mean with the
minimum and the maximum and no SD gives “mean_range”, where the mean is
kept and the SD is imputed from the range. A row with both a summary and
a sample mean and SD follows `prefer`, which is `"sample"` by default
because the sample statistics are the direct measurement. A row with
none of these, with a missing sample size, with values out of order, or
whose standard errors would be zero (a sample SD of zero, a summary
whose values are all equal) gets `NA` in every output column and is
named in one warning with its reason. The function does not stop, so the
whole table comes back. Output columns left by an earlier call are
replaced.

The factors are computed exactly at every sample size. At $`n = 4Q + 1`$
the reported quartiles are single order statistics. At any other sample
size each reported quantile is a fixed combination of two neighboring
order statistics, by R’s default definition (type 7 of
[`quantile()`](https://rdrr.io/r/stats/quantile.html)) unless the
argument `type` names another, and the factors there depend on that
definition. The closed-form shortcut formulas of Appendix A are within
half a percent of the exact factors at every sample size $`n = 4Q + 1`$
from 5 to 1001.

## Where to look afterwards

Each group’s output columns carry its suffix.

| Column | What it holds |
|----|----|
| `scenario` | `"S1"`, `"S2"`, `"S3"`, `"sample"` or `"mean_range"`, as read off the row |
| `mu_hat`, `sigma_hat` | the imputed mean and SD, or the sample ones |
| `se_naive` | $`\hat\sigma/\sqrt n`$, the standard error that treats the imputed mean as a sample mean |
| `se_mu` | $`C_\mu(n)\,\hat\sigma`$, the exact plug-in standard error of the imputed mean (for a sample mean, sd / sqrt(n)) |
| `se_sigma` | $`C_\sigma(n)\,\hat\sigma`$, the exact plug-in standard error of the imputed SD (for a sample SD, the standard deviation of a sample SD under normality with the SD plugged in) |
| `version_used` | `"blue"`, `"wan_luo_shi"` or `"sample"` |
| `exact` | `TRUE` for a row with an imputed value, whose factors are computed exactly at its sample size, `NA` for a sample row |

[`se_imputed_diff()`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_diff.md)
adds `yi`, the effect of the study, with `sei_naive` and `sei`, its
standard errors from the naive and from the exact standard errors, the
two groups’ variances added. By default the effect is the difference of
the two means (group 1 minus group 2), and `measure = "lnVR"` makes it
the log variability ratio, the log of the SD of group 1 over the SD of
group 2.
[`pool_imputed()`](https://zwulab.github.io/fivenumMeta/reference/pool_imputed.md)
pools these columns under both standard errors, and
`metafor::rma(yi, sei = sei)` takes them as they are.

The vignette `two-meta-analyses` reproduces the two meta-analyses of Wu,
Shi and Tong (2026): a meta-analysis of means in which studies report
the median and the range, on the data set `vitd` that ships with the
package, and a meta-analysis of SDs in which studies report the
quartiles, on this data set.

## Appendix A. The factors by scenario

[`inflation_factors()`](https://zwulab.github.io/fivenumMeta/reference/inflation_factors.md)
returns the ratio of the standard deviation of the imputed mean to that
of the sample mean, $`\kappa_\mu(n) = C_\mu(n)\sqrt n`$, and the ratio
of the standard deviation of the imputed SD to that of the sample SD,
$`\kappa_\sigma(n)`$, for the best linear unbiased estimator (BLUE, the
default) or the estimators of Wan, Luo and Shi
(`version = "wan_luo_shi"`).

``` r

grid <- c(5, 25, 101, 1001)
tab <- sapply(c("S1", "S2", "S3"),
              function(sc) sapply(grid, function(n) inflation_factors(n, sc)))
rownames(tab) <- paste0(rep(c("kappa_mu", "kappa_sigma"), length(grid)),
                        ", n = ", rep(grid, each = 2))
round(tab, 3)
#>                          S1    S2    S3
#> kappa_mu, n = 5       1.026 1.063 1.000
#> kappa_sigma, n = 5    1.089 1.683 1.070
#> kappa_mu, n = 25      1.131 1.075 1.034
#> kappa_sigma, n = 25   1.255 1.652 1.105
#> kappa_mu, n = 101     1.200 1.077 1.057
#> kappa_sigma, n = 101  1.705 1.650 1.253
#> kappa_mu, n = 1001    1.244 1.078 1.074
#> kappa_sigma, n = 1001 3.427 1.650 1.507
```

Under $`S_2`$ both factors are nearly constant in the sample size. Under
$`S_1`$ both grow with it: $`\kappa_\mu(n)`$ toward its limit
$`\sqrt{\pi/2}`$, about 25% above one, and $`\kappa_\sigma(n)`$ without
bound, because $`C_\sigma(n)`$ decays only like $`1/\log n`$ there while
the standard deviation of a sample SD decays like $`1/\sqrt n`$. Under
$`S_3`$ the extremes help at small sample sizes and the factors move
toward those of $`S_2`$ as the sample grows. The closed-form shortcut
formulas
[`c_mu_shortcut()`](https://zwulab.github.io/fivenumMeta/reference/c_mu_shortcut.md)
and
[`c_sigma_shortcut()`](https://zwulab.github.io/fivenumMeta/reference/c_sigma_shortcut.md)
give the factors $`C_\mu(n)`$ and $`C_\sigma(n)`$ without the
integration, for a reader without this package and for comparison:

``` r

c(exact = c_mu(25, "S2"), shortcut = c_mu_shortcut(25, "S2"))
#>     exact  shortcut 
#> 0.2150631 0.2149903
c(exact = c_sigma(25, "S2"), shortcut = c_sigma_shortcut(25, "S2"))
#>     exact  shortcut 
#> 0.2372150 0.2374137
```

## Appendix B. One study at a time

[`se_imputed()`](https://zwulab.github.io/fivenumMeta/reference/se_imputed.md)
takes one study’s reported values and sample size.

``` r

se_imputed(c(q1 = 10, med = 12, q3 = 15), n = 25, scenario = "S2")
#> Imputed mean and SD under S2, n = 25, the BLUE, quantile type 7
#>    mu_hat sigma_hat  se_naive     se_mu  se_sigma 
#>   12.3575    3.9252    0.7850    0.8442    0.9311 
#> Factors: C_mu(n) = 0.2151, C_sigma(n) = 0.2372. The weights are $w_mean and $w_sd; every field with str() or unclass().
```

The imputed mean is the median pulled a little toward the mid-quartile
range, the imputed SD is a multiple of the interquartile range, and the
exact standard error of the mean exceeds the naive one by the factor of
the table above at $`n = 25`$ under $`S_2`$.

## References

Balakrishnan N, Rychtar J, Taylor D, Walter SD (2022). Unified approach
to optimal estimation of mean and standard deviation from sample
summaries. Statistical Methods in Medical Research 31(11), 2087-2103.

Katzenschlager S, Zimmer AJ, Gottschalk C, Grafeneder J, Schmitz S,
Kraker S, Ganslmeier M, Muth A, Seitel A, Maier-Hein L, Benedetti A,
Larmann J, Weigand MA, McGrath S, Denkinger CM (2021). Can we predict
the severe course of COVID-19, a systematic review and meta-analysis of
indicators of clinical outcome? PLOS One 16, e0255154.

Luo D, Wan X, Liu J, Tong T (2018). Optimally estimating the sample mean
from the sample size, median, mid-range, and/or mid-quartile range.
Statistical Methods in Medical Research 27(6), 1785-1805.

McGrath S, Katzenschlager S, Zimmer AJ, Seitel A, Steele R, Benedetti A
(2023). Standard error estimation in meta-analysis of studies reporting
medians. Statistical Methods in Medical Research 32(2), 373-388.

McGrath S, Zhao X, Ozturk O, Katzenschlager S, Steele R, Benedetti A
(2024). metamedian: An R package for meta-analyzing studies reporting
medians. Research Synthesis Methods 15(2), 332-346.

Shi J, Luo D, Weng H, Zeng XT, Lin L, Chu H, Tong T (2020). Optimally
estimating the sample standard deviation from the five-number summary.
Research Synthesis Methods 11(5), 641-654.

Viechtbauer W (2010). Conducting meta-analyses in R with the metafor
package. Journal of Statistical Software 36(3), 1-48.

Wan X, Wang W, Liu J, Tong T (2014). Estimating the sample mean and
standard deviation from the sample size, median, range and/or
interquartile range. BMC Medical Research Methodology 14, 135.

Wu Z (2026). osmoments: Exact Moments of Order Statistics. R package
version 0.1.0.

Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard
deviation from a five-number summary. Manuscript.

Wu Z, Tong T (2026). osmoments: Exact means, variances and covariances
of order statistics in R. Manuscript.
