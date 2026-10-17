# Two meta-analysis examples: a mean case and a standard deviation (SD) case

``` r

library(fivenumMeta)
has_mm <- requireNamespace("metamedian", quietly = TRUE)  # the SD case's data
has_mf <- requireNamespace("metafor", quietly = TRUE)     # confirms the pooling
```

## What the two cases show

This vignette assumes that the reader knows what
[`se_imputed_diff()`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_diff.md)
and
[`se_imputed_df()`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_df.md)
return: each study’s imputed mean and SD with their naive and exact
standard errors, in the columns that
[`vignette("fivenumMeta")`](https://zwulab.github.io/fivenumMeta/articles/fivenumMeta.md)
describes. Start there if the columns, the scenarios or the conventions
are new.

Wu, Shi and Tong (2026) illustrate the exact standard errors on two
published meta-analyses, and this vignette runs both with the package.
In a meta-analysis of means in which three of six studies report the
median and the range, the exact standard errors raise the pooled
standard error by 7.1%. In a meta-analysis of SDs in which 29 of 51
studies report the quartiles, they raise the pooled standard error by
35% and bring $`I^2`$ from 86% to 70%. Both data sets are public. The
vitamin D studies ship with this package as the data set `vitd`, and the
age studies come with the `metamedian` package as its data set
`dat.age`. Every number in this vignette is computed when it is built.

Both cases pool with the common-effect model. Study $`i`$ gives an
effect estimate $`y_i`$ with variance $`v_i`$, and the model weights it
by $`w_i = 1/v_i`$:
``` math
\hat\theta = \frac{\sum_i w_i y_i}{\sum_i w_i}, \qquad
\mathrm{SE}(\hat\theta) = \Big(\sum_i w_i\Big)^{-1/2}, \qquad
\mathcal Q = \sum_i w_i (y_i - \hat\theta)^2, \qquad
I^2 = \max\Big\{0, \frac{\mathcal Q - (K - 1)}{\mathcal Q}\Big\},
```
where $`K`$ is the number of studies. $`\mathcal Q`$ and $`I^2`$ are the
heterogeneity measures of Higgins and Thompson (2002).

We run each meta-analysis twice: with the naive standard errors, which
treat the imputed mean and SD as a sample mean and SD, and with the
exact standard errors computed from the exact standard deviation
factors, which this package takes from the moments of the normal order
statistics of the `osmoments` package (Wu and Tong, 2026). That is, the
two analyses share the $`y_i`$ and differ only in the $`v_i`$. The
function
[`pool_imputed()`](https://zwulab.github.io/fivenumMeta/reference/pool_imputed.md)
computes the four quantities above twice from the output of
[`se_imputed_diff()`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_diff.md),
once with the naive and once with the exact standard errors, and prints
the two fits side by side. Both cases call it.

## A meta-analysis of means from medians and ranges

The data are the meta-analysis of serum vitamin D (nmol/L) in
tuberculosis cases against controls of Nnoaham and Clarke (2008), typed
from their Table 1 and shipped with this package as the data set `vitd`.
Three studies report the median, minimum and maximum of each group, two
a sample mean and SD, and one a mean with the minimum and maximum.

``` r

vitd
#>               study n.cases min.cases med.cases max.cases mean.cases sd.cases
#> 1        Davies, UK      40      2.25     16.00     74.25         NA       NA
#> 2 Grange, Indonesia      40     43.75     65.75    130.50         NA       NA
#> 3     Davies, Kenya      15     16.75     39.75     89.25         NA       NA
#> 4  Davies, Thailand      51        NA        NA        NA      69.50     24.5
#> 5   Chan, Hong Kong      22        NA        NA        NA      46.50     18.5
#> 6 Sasidharan, India      35      2.50        NA     75.00      26.75       NA
#>   n.controls min.controls med.controls max.controls mean.controls sd.controls
#> 1         40         9.00        27.25       132.50            NA          NA
#> 2         38        48.50        69.50       125.00            NA          NA
#> 3         15        26.25        65.50       114.75            NA          NA
#> 4         51           NA           NA           NA         95.50       29.25
#> 5         23           NA           NA           NA         52.25       15.75
#> 6         16        22.50           NA       145.00         48.50          NA
```

The study that reports a mean with the range (Sasidharan, India) has no
SD.
[`se_imputed_df()`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_df.md)
reads such a row as the scenario `mean_range`: the mean is kept as the
sample mean, and the SD is imputed from the two extremes alone, as the
range divided by the expected range of a normal sample of that size,
which is the best linear unbiased estimator (BLUE) of the SD from the
extremes. Its exact standard error is the standard deviation of that
estimator, from the moments of the two extremes. The standard error of
this study’s mean is then $`\hat\sigma/\sqrt n`$ under both
calculations, with the imputed SD plugged in.

The call takes the controls first, so the difference `yi` is controls
minus cases, as the paper reports it. Each group’s factors are taken at
its own sample size, exact at the even sizes 40 and 38, where the
reported median is the average of the two middle values under R’s
default quantile definition. The table lists each study’s scenario, its
difference with the two standard errors, and its weight in percent under
each.

``` r

d <- se_imputed_diff(vitd, suffixes = c(".controls", ".cases"))
fit_vd <- pool_imputed(d)
tab <- data.frame(study = d$study, scenario = d$scenario.cases,
                  yi = d$yi, se_naive = d$sei_naive, se_exact = d$sei,
                  w_naive = fit_vd$weights$naive,
                  w_exact = fit_vd$weights$exact)
print(tab, digits = 3, row.names = FALSE)
#>              study   scenario    yi se_naive se_exact w_naive w_exact
#>         Davies, UK         S1 15.42     5.23     6.01   19.30   16.77
#>  Grange, Indonesia         S1  3.04     4.30     4.93   28.57   24.89
#>      Davies, Kenya         S1 22.91     8.51     9.34    7.29    6.94
#>   Davies, Thailand     sample 26.00     5.34     5.34   18.49   21.20
#>    Chan, Hong Kong     sample  5.75     5.13     5.13   20.04   22.97
#>  Sasidharan, India mean_range 21.75     9.15     9.15    6.31    7.23
fit_vd
#> Pooled difference of means (group 1 minus group 2), common-effect model (EE), 6 studies
#>          estimate    se ci.lb ci.ub      pval     Q    I2
#> naive SE    12.85 2.298 8.346 17.35 2.240e-08 15.77 68.29
#> exact SE    13.34 2.460 8.517 18.16 5.876e-08 14.18 64.75
#> 95% confidence intervals and p-values from the normal distribution. Q has 5 degrees of freedom, and I2 is in percent.
#> The weights in percent are $weights. Every field with str() or unclass().
```

Compared with the naive calculation, the exact calculation raises the
standard errors of the three studies that report the median and the
range by 9.8% to 14.8%, the inflation of the imputed mean under that
scenario at these sample sizes. Their share of the weight falls from 55%
to 49%. The pooled difference moves from 12.8 to 13.3 nmol/L, and its
standard error from 2.30 to 2.46, a rise of 7.1%, so its confidence
interval widens by the same proportion.

## A meta-analysis of SDs from quartiles

A meta-analysis of variability studies whether an outcome is more
dispersed in one group than in another. Here the effect size is the log
variability ratio of Nakagawa et al. (2015), the log of the ratio of the
two groups’ SDs, $`y_i = \log(\hat\sigma_{i1}/\hat\sigma_{i2})`$. The
data are the 51 studies of the age of COVID-19 nonsurvivors (group 1)
against survivors (group 2) of Katzenschlager et al. (2021), shipped
with the `metamedian` package as its data set `dat.age`: 29 report the
three quartiles and 22 a sample mean and SD in both groups. The question
is whether age is less dispersed among those who died.
`measure = "lnVR"` gives $`y_i`$ with its two standard errors in place
of the difference of means.

``` r

d <- se_imputed_diff(metamedian::dat.age, measure = "lnVR")  # suffixes .g1, .g2
table(d$scenario.g1)
#> 
#>     S2 sample 
#>     29     22
fit_vr <- pool_imputed(d)
fit_vr
#> Pooled log variability ratio (SD of group 1 over SD of group 2), common-effect model (EE), 51 studies
#>          estimate      se   ci.lb   ci.ub      pval     Q    I2
#> naive SE  -0.2253 0.01610 -0.2568 -0.1937 1.768e-44 346.6 85.58
#> exact SE  -0.2152 0.02167 -0.2577 -0.1728 2.939e-23 165.3 69.75
#> 95% confidence intervals and p-values from the normal distribution. Q has 50 degrees of freedom, and I2 is in percent.
#> exp() of the estimate and of the limits gives the ratio of SDs.
#> The weights in percent are $weights. Every field with str() or unclass().
```

The variance of $`y_i`$ is approximated by the sum of one term per
group, to first order in each group’s SD (the delta method,
$`\mathrm{Var}(\log X) \approx \mathrm{Var}(X)/(\mathrm{E}X)^2`$). The
naive calculation gives every group the large-sample term of a sample SD
of $`n`$ normal observations, $`1/\{2(n - 1)\}`$. The exact calculation
keeps that term for a group that gives its sample SD and uses
$`C^2_\sigma(n)`$, the exact variance of the imputed SD in units of
$`\sigma^2`$, for a group whose SD is imputed from the quartiles, which
the output carries as `(se_sigma / sigma_hat)^2`:
``` math
v_i^{\mathrm{naive}} = \frac{1}{2(n_{i1} - 1)} + \frac{1}{2(n_{i2} - 1)}, \qquad
v_i^{\mathrm{exact}} = t(n_{i1}) + t(n_{i2}), \qquad
t(n) = \begin{cases}
  1/\{2(n - 1)\} & \text{for a sample SD,} \\
  C^2_\sigma(n) & \text{for an imputed SD.}
\end{cases}
```
The columns `sei_naive` and `sei` are the square roots of these two
variances. Both are first-order approximations of the variance of a log,
so the intervals below are normal approximations. The bias correction of
the log ratio of Nakagawa et al. (2015) is not applied, as in Wu, Shi
and Tong (2026). Both analyses use the plain log ratio.

For a group whose SD is imputed from the quartiles, $`C^2_\sigma(n)`$ is
2.34 to 2.72 times the naive term $`1/\{2(n - 1)\}`$. Under the exact
variance, the share of the weight on the 29 quartile studies falls from
72% to 49%. The pooled log ratio changes little, from -0.225 to -0.215,
so the SD of age among nonsurvivors is about 0.80 times that among
survivors either way. Its standard error, however, rises by 35%, from
0.0161 to 0.0217, and its confidence interval widens by the same
proportion. $`\mathcal Q`$ falls from 346.6 to 165.3, and $`I^2`$ from
86% to 70%. Since the two fits share the $`y_i`$ and differ only in the
variances, the fall in $`\mathcal Q`$ and $`I^2`$ comes from the larger
within-study variances alone: the extra spread of an imputed SD, which
the naive variance leaves out, was counted as heterogeneity.

With the `metafor` package, the exact fit is
`rma(yi, sei = sei, data = d, method = "EE")`, which returns the same
pooled estimate, standard error and $`\mathcal Q`$, and the naive fit
takes `sei_naive` in the same way:

``` r

f <- metafor::rma(yi, sei = sei, data = d, method = "EE")
c(estimate = as.numeric(f$b), se = f$se, Q = f$QE, I2 = f$I2)
#>     estimate           se            Q           I2 
#>  -0.21524814   0.02166619 165.27554705  69.74749085
```

## References

Higgins JPT, Thompson SG (2002). Quantifying heterogeneity in a
meta-analysis. Statistics in Medicine 21(11), 1539-1558.

Katzenschlager S, Zimmer AJ, Gottschalk C, Grafeneder J, Schmitz S,
Kraker S, Ganslmeier M, Muth A, Seitel A, Maier-Hein L, Benedetti A,
Larmann J, Weigand MA, McGrath S, Denkinger CM (2021). Can we predict
the severe course of COVID-19, a systematic review and meta-analysis of
indicators of clinical outcome? PLOS One 16, e0255154.

McGrath S, Zhao X, Ozturk O, Katzenschlager S, Steele R, Benedetti A
(2024). metamedian: An R package for meta-analyzing studies reporting
medians. Research Synthesis Methods 15(2), 332-346.

Nakagawa S, Poulin R, Mengersen K, Reinhold K, Engqvist L, Lagisz M,
Senior AM (2015). Meta-analysis of variation: ecological and
evolutionary applications and beyond. Methods in Ecology and Evolution
6(2), 143-152.

Nnoaham KE, Clarke A (2008). Low serum vitamin D levels and
tuberculosis: a systematic review and meta-analysis. International
Journal of Epidemiology 37(1), 113-119.

Viechtbauer W (2010). Conducting meta-analyses in R with the metafor
package. Journal of Statistical Software 36(3), 1-48.

Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard
deviation from a five-number summary. Manuscript.

Wu Z, Tong T (2026). osmoments: Exact means, variances and covariances
of order statistics in R. Manuscript.
