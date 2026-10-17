# Serum vitamin D in tuberculosis cases and controls: six studies

The six studies of the meta-analysis of serum vitamin D (nmol/L) in
tuberculosis cases against controls of Nnoaham and Clarke (2008), typed
from their Table 1. Three studies report the median, minimum and maximum
of each group (Davies, UK; Grange, Indonesia; Davies, Kenya: the
scenario `"S1"` of
[`se_imputed_df`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_df.md)),
two a sample mean and SD (Davies, Thailand; Chan, Hong Kong:
`"sample"`), and one a mean with the minimum and maximum (Sasidharan,
India: `"mean_range"`, the mean kept and the SD imputed from the range).
The same six studies are the real-data example of Luo et al. (2018) and
of Yang, Hutson and Wang (2022). The vignette
[`vignette("two-meta-analyses")`](https://zwulab.github.io/fivenumMeta/articles/two-meta-analyses.md)
pools them under the naive and the exact standard errors, as Wu, Shi and
Tong (2026) do.

## Usage

``` r
vitd
```

## Format

A data frame with 6 rows and 13 columns: `study`, the first author and
the country, as the source labels the study, then for each group, with
the suffix `.cases` or `.controls`, the columns `n` (the sample size),
`min`, `med`, `max` (the minimum, the median, the maximum), `mean` and
`sd` (the sample mean and SD), `NA` where the study does not report the
value.

## Source

Nnoaham KE, Clarke A (2008). Low serum vitamin D levels and
tuberculosis: a systematic review and meta-analysis. International
Journal of Epidemiology 37(1), 113-119, Table 1. Chan et al. has 22
cases and 23 controls there, in Table 5 of Yang et al. (2022) and in
Tables 2 and 3 of Luo et al. (2018), whose Table 1 prints 24 and 24.

## References

Luo D, Wan X, Liu J, Tong T (2018). Optimally estimating the sample mean
from the sample size, median, mid-range, and/or mid-quartile range.
Statistical Methods in Medical Research 27(6), 1785-1805.

Yang X, Hutson AD, Wang D (2022). A generalized BLUE approach for
combining location and scale information in a meta-analysis. Journal of
Applied Statistics 49(15), 3846-3867.

Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard
deviation from a five-number summary. Manuscript.

## Examples

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
# the difference of means, controls minus cases, with its naive and its exact
# standard error, and the pooled difference under each
d <- se_imputed_diff(vitd, suffixes = c(".controls", ".cases"))
d[, c("study", "scenario.cases", "yi", "sei_naive", "sei")]
#>               study scenario.cases        yi sei_naive      sei
#> 1        Davies, UK             S1 15.424644  5.230349 6.006378
#> 2 Grange, Indonesia             S1  3.039868  4.298886 4.930551
#> 3     Davies, Kenya             S1 22.912639  8.508270 9.339537
#> 4  Davies, Thailand         sample 26.000000  5.342784 5.342784
#> 5   Chan, Hong Kong         sample  5.750000  5.132460 5.132460
#> 6 Sasidharan, India     mean_range 21.750000  9.145621 9.145621
pool_imputed(d)
#> Pooled difference of means (group 1 minus group 2), common-effect model (EE), 6 studies
#>          estimate    se ci.lb ci.ub      pval     Q    I2
#> naive SE    12.85 2.298 8.346 17.35 2.240e-08 15.77 68.29
#> exact SE    13.34 2.460 8.517 18.16 5.876e-08 14.18 64.75
#> 95% confidence intervals and p-values from the normal distribution. Q has 5 degrees of freedom, and I2 is in percent.
#> The weights in percent are $weights. Every field with str() or unclass().
```
