# The effect of two groups per study, with its naive and exact standard errors

For a data set with two groups side by side (the columns of one group
carry one suffix, those of the other group another), runs
[`se_imputed_df`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_df.md)
on each group and appends the effect of the study, group 1 against group
2, with its standard error under the naive and under the exact standard
errors: the difference of means (`measure = "MD"`, the default) or the
log variability ratio (`measure = "lnVR"`), in the column names that
[`metafor::rma()`](https://wviechtb.github.io/metafor/reference/rma.uni.html)
takes as they are, `yi` and `sei`.
[`pool_imputed`](https://zwulab.github.io/fivenumMeta/reference/pool_imputed.md)
pools the result in one call.

## Usage

``` r
se_imputed_diff(
  data,
  suffixes = c(".g1", ".g2"),
  ...,
  measure = c("MD", "lnVR"),
  version = c("blue", "wan_luo_shi"),
  prefer = c("sample", "imputed"),
  type = 7
)
```

## Arguments

- data:

  A data frame, one row per study.

- suffixes:

  The two suffixes of the groups' columns, group 1 first. The default is
  the layout of the data sets shipped with the metamedian package,
  `".g1"` and `".g2"`.

- ...:

  The column-name arguments of
  [`se_imputed_df`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_df.md)
  (`n`, `min`, `q1`, `med`, `q3`, `max`, `mean`, `sd`), passed to both
  groups.

- measure:

  The effect per study: `"MD"` (default), the difference of the two
  means, group 1 minus group 2; or `"lnVR"`, the log variability ratio
  of Nakagawa et al. (2015), the log of the SD of group 1 over the SD of
  group 2. See Details.

- version:

  `"blue"` (default) or `"wan_luo_shi"`, as in
  [`se_imputed`](https://zwulab.github.io/fivenumMeta/reference/se_imputed.md).

- prefer:

  What to use for a row that has both a summary and a sample mean and
  SD: `"sample"` (default), the sample statistics, or `"imputed"`.

- type:

  The definition of the reported sample quantiles, 1 to 9 as in
  [`quantile`](https://rdrr.io/r/stats/quantile.html); 7, R's default.
  One value for every row.

## Value

`data` with the output columns of
[`se_imputed_df`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_df.md)
for each group (each with its suffix) and three more: `yi`, the effect;
`sei_naive`, its standard error from the naive standard errors; and
`sei`, its standard error from the exact ones. The result carries the
attribute `"measure"`, which
[`pool_imputed`](https://zwulab.github.io/fivenumMeta/reference/pool_imputed.md)
reads for its heading. A row that
[`se_imputed_df()`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_df.md)
could not use in either group has `NA` in the three columns.

## Details

The two groups are independent samples, so the variances add.

Under `"MD"`, \\y_i = \hat\mu\_{i1} - \hat\mu\_{i2}\\, `sei_naive` \\=
(\mathrm{se\\naive}\_{i1}^2 + \mathrm{se\\naive}\_{i2}^2)^{1/2}\\ and
`sei` \\= (\mathrm{se\\mu}\_{i1}^2 + \mathrm{se\\mu}\_{i2}^2)^{1/2}\\.

Under `"lnVR"`, \\y_i = \log(\hat\sigma\_{i1}/\hat\sigma\_{i2})\\, and
the variance of \\y_i\\ is approximated by the sum of one term per
group, to first order in each group's SD (the delta method,
\\\mathrm{Var}(\log X) \approx \mathrm{Var}(X)/(\mathrm{E}X)^2\\). The
naive variance gives every group the large-sample term of the log of a
sample SD of \\n\\ normal observations, \\1/\\2(n - 1)\\\\, so
`sei_naive` \\= \\1/(2(n\_{i1} - 1)) + 1/(2(n\_{i2} - 1))\\^{1/2}\\. The
exact variance keeps that term for a group that gives its sample SD and
uses \\C^2\_\sigma(n)\\, the exact variance of the imputed SD in units
of \\\sigma^2\\, for a group whose SD is imputed (from the quartiles,
the extremes, or the range alone), which the group's columns carry as
`(se_sigma / sigma_hat)^2`; `sei` is the square root of the sum of the
two terms. These are the two variances of Wu, Shi and Tong (2026): exact
in the factor \\C\_\sigma(n)\\, first-order in the log, so neither is
the finite-sample variance of the log ratio, and the intervals of
[`pool_imputed`](https://zwulab.github.io/fivenumMeta/reference/pool_imputed.md)
are normal approximations. The small-sample bias correction of Nakagawa
et al. (2015), which `metafor::escalc(measure = "VR")` adds to \\y_i\\,
is not applied, as in Wu, Shi and Tong (2026). Both analyses use the
plain log ratio.

## References

Nakagawa S, Poulin R, Mengersen K, Reinhold K, Engqvist L, Lagisz M,
Senior AM (2015). Meta-analysis of variation: ecological and
evolutionary applications and beyond. Methods in Ecology and Evolution
6(2), 143-152.

Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard
deviation from a five-number summary. Manuscript.

## See also

[`se_imputed_df`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_df.md),
[`se_imputed`](https://zwulab.github.io/fivenumMeta/reference/se_imputed.md),
[`pool_imputed`](https://zwulab.github.io/fivenumMeta/reference/pool_imputed.md).

## Examples

``` r
two <- data.frame(n.g1 = c(25, 30), q1.g1 = c(10, 11), med.g1 = c(12, 13), q3.g1 = c(15, 16),
                  n.g2 = c(29, 40), q1.g2 = c(8, 9), med.g2 = c(10, 11), q3.g2 = c(13, 14))
d <- se_imputed_diff(two)
d[, c("yi", "sei_naive", "sei")]
#>         yi sei_naive       sei
#> 1 2.001110 1.0672925 1.1478665
#> 2 2.002044 0.9333307 0.9977329
pool_imputed(d)                                    # the pooled difference, both SEs
#> Pooled difference of means (group 1 minus group 2), common-effect model (EE), 2 studies
#>          estimate     se  ci.lb ci.ub     pval         Q I2
#> naive SE    2.002 0.7026 0.6246 3.379 0.004386 4.337e-07  0
#> exact SE    2.002 0.7530 0.5257 3.478 0.007858 3.769e-07  0
#> 95% confidence intervals and p-values from the normal distribution. Q has 1 degrees of freedom, and I2 is in percent.
#> The weights in percent are $weights. Every field with str() or unclass().
v <- se_imputed_diff(two, measure = "lnVR")       # the log ratio of the SDs
v[, c("yi", "sei_naive", "sei")]
#>            yi sei_naive       sei
#> 1 0.008045533 0.1966989 0.3233439
#> 2 0.011700746 0.1733837 0.2791642
# columns named otherwise, two groups with the suffixes "_t" and "_c"
spread <- data.frame(N_t = 25, Min_t = 3, Median_t = 12, Max_t = 30,
                     N_c = 29, Min_c = 2, Median_c = 10, Max_c = 25)
se_imputed_diff(spread, suffixes = c("_t", "_c"), n = "N", min = "Min",
                med = "Median", max = "Max")[, c("yi", "sei_naive", "sei")]
#>         yi sei_naive      sei
#> 1 2.337257  1.730794 1.963484
if (FALSE) { # \dontrun{
metafor::rma(yi, sei = sei, data = d)          # the pooled difference with the exact SEs
metafor::rma(yi, sei = sei_naive, data = d)    # the same with the naive SEs, for comparison
} # }
```
