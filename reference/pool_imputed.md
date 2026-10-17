# Pool the studies under the naive and the exact standard errors, side by side

Takes the output of
[`se_imputed_diff`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_diff.md)
(or any data frame with the columns `yi`, `sei_naive` and `sei`) and
pools the studies twice, once with the naive standard errors and once
with the exact ones, by the common-effect model (`method = "EE"`, the
default) or the DerSimonian-Laird random-effects model
(`method = "DL"`). The two fits share the effect estimates and differ in
the variances alone, so their difference shows what the exact standard
errors change: the pooled estimate, its standard error and confidence
interval, the heterogeneity statistic \\\mathcal Q\\ and the index
\\I^2\\, and the weights of the studies.

## Usage

``` r
pool_imputed(data, method = c("EE", "DL"), level = 0.95)

# S3 method for class 'pool_imputed'
print(x, digits = 4L, ...)
```

## Arguments

- data:

  A data frame with the columns `yi` (the effect of each study),
  `sei_naive` (its standard error from the naive standard errors) and
  `sei` (from the exact ones), as
  [`se_imputed_diff`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_diff.md)
  returns it; the three must be numeric. A row whose effect is not
  finite, or whose standard error in either column is not finite and
  positive, is dropped from both fits, so that the two rows compare the
  same studies, with a warning that names the row and the column.

- method:

  `"EE"` (default), the common-effect (equal-effects) model, or `"DL"`,
  the random-effects model with the DerSimonian-Laird estimator of
  \\\tau^2\\. The codes are those of
  [`metafor::rma()`](https://wviechtb.github.io/metafor/reference/rma.uni.html).

- level:

  The confidence level of the interval, a number strictly between 0 and
  1; 0.95 by default.

- x:

  An object returned by `pool_imputed()`.

- digits:

  The number of significant digits the print shows.

- ...:

  Passed to
  [`print.data.frame`](https://rdrr.io/r/base/print.dataframe.html) for
  the table.

## Value

An object of class `"pool_imputed"`: a list, used as any list, whose
print shows the heading and the two-row table. Its components: `fits`, a
data frame with the rows `naive` and `exact` and the columns `estimate`,
`se`, `ci.lb`, `ci.ub`, `pval`, `Q`, `I2` (in percent) and `tau2`;
`weights`, a data frame with one row per study used and the columns
`naive` and `exact`, the weights in percent of their total; `k`, the
number of studies used; `method`; `level`; and `measure`, read from the
attribute that
[`se_imputed_diff`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_diff.md)
sets (`"MD"` or `"lnVR"`), `NA` when `data` carries none (a column
subset `data[, cols]` drops the attribute, a row subset keeps it;
without it the print says "Pooled estimate").

## Details

With \\v_i\\ the squared standard error of study \\i\\, \\w_i = 1/v_i\\
and \\K\\ studies, the common-effect model gives \$\$\hat\theta =
\frac{\sum_i w_i y_i}{\sum_i w_i}, \qquad \mathrm{SE}(\hat\theta) =
\Big(\sum_i w_i\Big)^{-1/2},\$\$ \$\$\mathcal Q = \sum_i w_i (y_i -
\hat\theta)^2, \qquad I^2 = \max\Big\\0, \frac{\mathcal Q - (K -
1)}{\mathcal Q}\Big\\,\$\$ with \\\tau^2 = 0\\. \\\mathcal Q\\ and
\\I^2\\ are the heterogeneity measures of Higgins and Thompson (2002):
when the studies share one effect and the \\v_i\\ are the variances of
the \\y_i\\, \\\mathcal Q\\ has expectation \\K - 1\\, so \\I^2\\ is the
share of \\\mathcal Q\\ beyond that. The DerSimonian-Laird model
(DerSimonian and Laird, 1986) keeps \\\mathcal Q\\ and \\I^2\\ from
these weights and sets \$\$\hat\tau^2 = \max\Big\\0, \frac{\mathcal Q -
(K - 1)}{\sum_i w_i - \sum_i w_i^2 / \sum_i w_i}\Big\\,\$\$ \$\$w_i^\ast
= \frac{1}{v_i + \hat\tau^2},\$\$ and takes \\\hat\theta\\ and
\\\mathrm{SE}(\hat\theta)\\ as above with \\w_i^\ast\\ in place of
\\w_i\\. The confidence interval is \\\hat\theta \pm z\_{1 -
\alpha/2}\\\mathrm{SE}(\hat\theta)\\ and the p-value is the two-sided
normal test of \\\theta = 0\\, as the default of
[`metafor::rma()`](https://wviechtb.github.io/metafor/reference/rma.uni.html).
The weights are \\w_i\\ under `"EE"` and \\w_i^\ast\\ under `"DL"`, in
percent of their total.

Everything is computed twice, with \\v_i\\ the square of `sei_naive` and
with \\v_i\\ the square of `sei`. Under the naive standard errors an
imputed mean or SD is treated as a sample mean or SD, so the variances
of the studies that report a summary are too small; the exact row shows
the fit with the variances of Wu, Shi and Tong (2026). Under the
common-effect model the pooled standard error is \\(\sum_i w_i)^{-1/2}\\
and \\\mathcal Q\\ is the minimum over \\\theta\\ of \\\sum_i w_i (y_i -
\theta)^2\\, so when no exact variance is below its naive one (the case
for the difference of means, since \\\kappa\_\mu(n) \ge 1\\) the exact
row has the larger pooled standard error and the smaller \\\mathcal Q\\
and \\I^2\\. Under `"DL"` the direction of the pooled standard error is
not fixed, because the re-estimated \\\hat\tau^2\\ enters every weight
and a smaller \\\hat\tau^2\\ can outweigh the larger \\v_i\\.

What this function does not supply: a p-value for \\\mathcal Q\\, REML
and the other estimators of \\\tau^2\\, the Knapp-Hartung adjustment and
moderators, which `metafor::rma(yi, sei = sei, data = d)` fits from the
same columns (Viechtbauer, 2010). The exact reference distribution of
the pooled estimate under imputed standard errors is the subject of a
later version: the intervals and p-values here are normal
approximations, as metafor's default ones are.

## References

DerSimonian R, Laird N (1986). Meta-analysis in clinical trials.
Controlled Clinical Trials 7(3), 177-188.

Higgins JPT, Thompson SG (2002). Quantifying heterogeneity in a
meta-analysis. Statistics in Medicine 21(11), 1539-1558.

Viechtbauer W (2010). Conducting meta-analyses in R with the metafor
package. Journal of Statistical Software 36(3), 1-48.

Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard
deviation from a five-number summary. Manuscript.

## See also

[`se_imputed_diff`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_diff.md),
which builds `data`;
[`se_imputed_df`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_df.md).

## Examples

``` r
# a meta-analysis of means: controls minus cases, the data set vitd
d <- se_imputed_diff(vitd, suffixes = c(".controls", ".cases"))
fit <- pool_imputed(d)
fit
#> Pooled difference of means (group 1 minus group 2), common-effect model (EE), 6 studies
#>          estimate    se ci.lb ci.ub      pval     Q    I2
#> naive SE    12.85 2.298 8.346 17.35 2.240e-08 15.77 68.29
#> exact SE    13.34 2.460 8.517 18.16 5.876e-08 14.18 64.75
#> 95% confidence intervals and p-values from the normal distribution. Q has 5 degrees of freedom, and I2 is in percent.
#> The weights in percent are $weights. Every field with str() or unclass().
fit$fits["exact", c("estimate", "se")]
#>       estimate       se
#> exact 13.33852 2.459836
round(fit$weights, 1)                      # each study's weight, percent
#>   naive exact
#> 1  19.3  16.8
#> 2  28.6  24.9
#> 3   7.3   6.9
#> 4  18.5  21.2
#> 5  20.0  23.0
#> 6   6.3   7.2
pool_imputed(d, method = "DL", level = 0.9)
#> Pooled difference of means (group 1 minus group 2), random-effects model, DerSimonian-Laird (DL), 6 studies
#>          estimate    se ci.lb ci.ub      pval     Q    I2  tau2
#> naive SE    14.73 4.266 7.713 21.75 0.0005542 15.77 68.29 71.27
#> exact SE    14.75 4.312 7.660 21.84 0.0006230 14.18 64.75 69.26
#> 90% confidence intervals and p-values from the normal distribution. Q has 5 degrees of freedom, and I2 is in percent.
#> The weights in percent are $weights. Every field with str() or unclass().
if (FALSE) { # \dontrun{
# a meta-analysis of SDs: the log variability ratio on metamedian::dat.age
a <- se_imputed_diff(metamedian::dat.age, measure = "lnVR")
pool_imputed(a)
metafor::rma(yi, sei = sei, data = a, method = "EE")   # the same exact fit
} # }
```
