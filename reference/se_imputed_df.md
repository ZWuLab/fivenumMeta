# Imputed means and SDs with their standard errors for every study of a data frame

The data-frame interface of
[`se_imputed`](https://zwulab.github.io/fivenumMeta/reference/se_imputed.md)
for a meta-analysis: one row per study, with the columns that such a
data set usually carries (the sample size, the minimum, the first
quartile, the median, the third quartile, the maximum, and the sample
mean and SD where a study gives them). For each row the function reads
which values are present, decides the data-reporting scenario, imputes
the mean and the SD where they are not given, and appends them with
their standard errors: the naive \\\hat\sigma/\sqrt n\\ that treats the
imputed mean as a sample mean, and the plug-in standard errors
\\C\_\mu(n)\hat\sigma\\ and \\C\_\sigma(n)\hat\sigma\\ built from the
factors of Wu, Shi and Tong (2026).

## Usage

``` r
se_imputed_df(
  data,
  n = "n",
  min = "min",
  q1 = "q1",
  med = "med",
  q3 = "q3",
  max = "max",
  mean = "mean",
  sd = "sd",
  suffix = "",
  version = c("blue", "wan_luo_shi"),
  prefer = c("sample", "imputed"),
  append = TRUE,
  type = 7
)
```

## Arguments

- data:

  A data frame, one row per study.

- n, min, q1, med, q3, max, mean, sd:

  The names of the columns holding the sample size, the minimum, the
  first quartile, the median, the third quartile, the maximum, and the
  sample mean and SD, each with `suffix` appended. An argument set to
  `NULL`, or naming a column that `data` does not have, counts as absent
  for every row. A missing value (`NA`) counts as absent for that row.
  The columns must be numeric; a column that is all `NA` may be logical,
  as [`read.csv()`](https://rdrr.io/r/utils/read.table.html) reads such
  a column.

- suffix:

  A string appended to every column name above and to every output
  column, for data sets that carry two groups side by side (`".g1"`,
  `".g2"`). The default is the empty string.

- version:

  `"blue"` (default) or `"wan_luo_shi"`, as in
  [`se_imputed`](https://zwulab.github.io/fivenumMeta/reference/se_imputed.md).

- prefer:

  What to use for a row that has both a summary and a sample mean and
  SD: `"sample"` (default), the sample statistics, or `"imputed"`.

- append:

  `TRUE` (default) returns `data` with the output columns appended;
  `FALSE` returns the output columns alone.

- type:

  The definition of the reported sample quantiles, 1 to 9 as in
  [`quantile`](https://rdrr.io/r/stats/quantile.html); 7, R's default.
  One value for every row.

## Value

A data frame. The output columns, each with `suffix` appended, are
`scenario` (`"S1"`, `"S2"`, `"S3"`, `"sample"` or `"mean_range"`, see
Details), `mu_hat`, `sigma_hat`, `se_naive` (\\\hat\sigma/\sqrt n\\),
`se_mu` (\\C\_\mu(n)\hat\sigma\\ for an imputed mean, \\\hat\sigma/\sqrt
n\\ for a sample mean), `se_sigma` (\\C\_\sigma(n)\hat\sigma\\ for an
imputed SD, and for a sample SD \\\mathrm{sd}\sqrt{1 - b_S(n)^2}\\, the
standard deviation of a sample SD under normality with the SD plugged
in), `version_used` (`"blue"`, `"wan_luo_shi"`, or `"sample"` when
nothing is imputed) and `exact` (`TRUE` for a row with an imputed value,
whose factors are computed exactly at its \\n\\, `NA` for a `"sample"`
row).

## Details

The scenario of a row is read off its non-missing values, in this order.
A sample mean and SD give `"sample"`, the two used as they are, unless
`prefer = "imputed"` and the row also has a complete summary. All five
values give `"S3"`; the first quartile, the median and the third
quartile give `"S2"` (a single extreme beside them is ignored); the
minimum, the median and the maximum give `"S1"`. A mean with the minimum
and the maximum, and no SD, gives `"mean_range"` when no summary is
complete: the mean is kept as the sample mean, and the SD is imputed
from the two extremes alone, as the range divided by its expected value
in a normal sample of size \\n\\ (the BLUE of \\\sigma\\ from the
extremes, under `version = "blue"`) or by the closed-form divisor
\\\xi(n)\\ of the Wan-Luo-Shi version (under `"wan_luo_shi"`), with
`se_sigma` \\= C\_\sigma(n)\hat\sigma\\ from the exact moments of the
two extremes. These are the SD estimators of `"S1"` for a sample whose
median sits at the center of the order statistics, where the median
carries no weight, computed here without a median, so the quantile
definition `type` plays no part. Such a row has `se_naive` and `se_mu`
both equal to \\\hat\sigma/\sqrt n\\, the standard error of a sample
mean with the imputed SD plugged in, and needs \\n \ge 3\\. A mean
outside its range gets `NA`, as an error in the row. A mean without an
SD beside a complete summary is ignored, and the summary gives both
estimates.

A row with none of these, with a missing, infinite or non-integer sample
size, with a non-finite mean, with a non-finite, negative or zero SD,
with a summary whose values are all equal or a range of zero (its
imputed SD is zero), with values out of order or a mean outside its
reported range, gets `NA` in every output column and is named in one
warning with its reason; the function does not stop, so that the whole
table comes back. A standard error of zero would give the study an
infinite inverse-variance weight, which is why such a row is set to `NA`
rather than kept. Output columns that an earlier call left in `data` are
replaced, so a data frame can be run again after a value, the `version`
or the `type` changes.

Every imputed row is computed exactly at its own \\n\\, as in
[`se_imputed`](https://zwulab.github.io/fivenumMeta/reference/se_imputed.md):
at a sample size not of the form \\4Q + 1\\ the reported quartiles are
interpolated under definition `type`.

Assumptions: independent normal data within each study. The standard
errors are plug-in quantities, \\\sigma\\ replaced by \\\hat\sigma\\;
their reference distribution is not Student's t and is not supplied by
this version. A study whose data are visibly skewed (the quartiles far
from symmetric around the median) should be screened before the
normal-model estimators are applied, as the source papers of those
estimators recommend.

## References

Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard
deviation from a five-number summary. Manuscript.

## See also

[`se_imputed`](https://zwulab.github.io/fivenumMeta/reference/se_imputed.md)
for one study,
[`se_imputed_diff`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_diff.md)
for two groups per row,
[`pool_imputed`](https://zwulab.github.io/fivenumMeta/reference/pool_imputed.md)
to pool the studies.

## Examples

``` r
studies <- data.frame(n = c(25, 30, 41, 25, 35),
                      min = c(3, NA, 1, NA, 2), q1 = c(10, 10, 9, NA, NA),
                      med = c(12, 12, 12, NA, NA), q3 = c(15, 15, 16, NA, NA),
                      max = c(30, NA, 30, NA, 28),
                      mean = c(NA, NA, NA, 13, 12), sd = c(NA, NA, NA, 3, NA))
se_imputed_df(studies)   # S3, S2, S1, a sample mean and SD, a mean with the range
#>    n min q1 med q3 max mean sd   scenario   mu_hat sigma_hat  se_naive
#> 1 25   3 10  12 15  30   NA NA         S3 13.02497  5.907023 1.1814046
#> 2 30  NA 10  12 15  NA   NA NA         S2 12.35952  3.883653 0.7090548
#> 3 41   1  9  12 16  30   NA NA         S3 12.72772  6.168222 0.9633144
#> 4 25  NA NA  NA NA  NA   13  3     sample 13.00000  3.000000 0.6000000
#> 5 35   2 NA  NA NA  28   12 NA mean_range 12.00000  6.171054 1.0430985
#>       se_mu  se_sigma version_used exact
#> 1 1.2219709 0.9367854         blue  TRUE
#> 2 0.7573922 0.8191412         blue  TRUE
#> 3 1.0053976 0.7902086         blue  TRUE
#> 4 0.6000000 0.4307056       sample    NA
#> 5 1.0430985 0.9957990         blue  TRUE
se_imputed_df(studies, version = "wan_luo_shi", append = FALSE)
#>     scenario   mu_hat sigma_hat  se_naive     se_mu  se_sigma version_used
#> 1         S3 13.02860  5.914788 1.1829577 1.2235941 0.9387881  wan_luo_shi
#> 2         S2 12.35650  3.892101 0.7105972 0.7590482 0.8227088  wan_luo_shi
#> 3         S3 12.72172  6.175581 0.9644636 1.0066064 0.7924849  wan_luo_shi
#> 4     sample 13.00000  3.000000 0.6000000 0.6000000 0.4307056       sample
#> 5 mean_range 12.00000  6.181485 1.0448616 1.0448616 0.9991682  wan_luo_shi
#>   exact
#> 1  TRUE
#> 2  TRUE
#> 3  TRUE
#> 4    NA
#> 5  TRUE
# columns named as a spreadsheet often names them
names(studies) <- c("N", "Min", "Q1", "Median", "Q3", "Max", "Mean", "SD")
se_imputed_df(studies, n = "N", min = "Min", q1 = "Q1", med = "Median", q3 = "Q3",
              max = "Max", mean = "Mean", sd = "SD", append = FALSE)
#>     scenario   mu_hat sigma_hat  se_naive     se_mu  se_sigma version_used
#> 1         S3 13.02497  5.907023 1.1814046 1.2219709 0.9367854         blue
#> 2         S2 12.35952  3.883653 0.7090548 0.7573922 0.8191412         blue
#> 3         S3 12.72772  6.168222 0.9633144 1.0053976 0.7902086         blue
#> 4     sample 13.00000  3.000000 0.6000000 0.6000000 0.4307056       sample
#> 5 mean_range 12.00000  6.171054 1.0430985 1.0430985 0.9957990         blue
#>   exact
#> 1  TRUE
#> 2  TRUE
#> 3  TRUE
#> 4    NA
#> 5  TRUE
```
