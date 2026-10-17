# Inflation factors of the imputed mean and SD

The inflation factors (IF) of Wu, Shi and Tong (2026): the ratio of the
standard deviation of the imputed mean to that of the sample mean,
\\\kappa\_\mu(n) = \mathrm{SD}(\hat\mu)/\mathrm{SD}(\bar X) =
C\_\mu(n)\sqrt n\\, and of the imputed SD to that of the sample SD
\\S\\, \\\kappa\_\sigma(n) = \mathrm{SD}(\hat\sigma)/\mathrm{SD}(S) =
C\_\sigma(n)/\sqrt{1 - b_S(n)^2}\\, where \\b_S(n) =
\sqrt{2/(n-1)}\\\Gamma(n/2)/\Gamma\\(n-1)/2\\\\ is the bias factor of
\\S\\, \\E\[S\] = b_S(n)\sigma\\. The factors \\C\_\mu(n)\\ and
\\C\_\sigma(n)\\ are the exact ones of
[`c_mu`](https://zwulab.github.io/fivenumMeta/reference/c_mu.md) and
[`c_sigma`](https://zwulab.github.io/fivenumMeta/reference/c_sigma.md),
from one exact moment set.

## Usage

``` r
inflation_factors(
  n,
  scenario = c("S3", "S1", "S2", "S4", "S5"),
  version = c("blue", "wan_luo_shi"),
  type = 7
)
```

## Arguments

- n:

  Integer sample size: at least 3 under S1 and S2, 5 under S3 and S4,
  and 7 under S5, under the default quantile definition. Any such `n` is
  computed exactly; see Details. Another definition may need a larger
  \\n\\: where two reported values are the same order statistic (type 3
  at \\n = 5\\ under S3, for example) the BLUE is not defined and is
  refused.

- scenario:

  Character. One of `"S1"` (min, median, max), `"S2"` (three quartiles),
  `"S3"` (all five). McGrath 2023 labels; note Wan 2014's C2/C3 are
  S3/S2 respectively. Also `"S4"` (the quartiles with the 10th and 90th
  percentiles) and `"S5"` (all five with the two percentiles), the names
  of Wu, Shi and Tong (2026) for Scenarios 5 and 7 of Balakrishnan et
  al. (2022); only the BLUE (`"blue"`) is defined for these two.

- version:

  `"blue"` (default) for the BLUE of both, or `"wan_luo_shi"` for the
  estimators of Wan, Luo and Shi (Luo's mean; Wan's SD under S1 and S2,
  Shi's under S3). Under S4 and S5 only `"blue"` is defined.

- type:

  The definition of the sample quantiles, an integer from 1 to 9 as in
  [`quantile`](https://rdrr.io/r/stats/quantile.html); 7, R's default,
  unless a study states another.

## Value

Named numeric vector: `kappa_mu`, `kappa_sigma`.

## Details

Assumptions as in
[`imputed_moments`](https://zwulab.github.io/fivenumMeta/reference/imputed_moments.md):
iid normal data, the reported values following definition `type`. An IF
of 1.2 means a 20% larger SD than that of the full-sample statistic.

## References

Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard
deviation from a five-number summary. Manuscript.

## Examples

``` r
inflation_factors(41, "S2")                       # the BLUE: about 1.076 and 1.651
#>    kappa_mu kappa_sigma 
#>    1.076399    1.651035 
inflation_factors(41, "S2", version = "wan_luo_shi")
#>    kappa_mu kappa_sigma 
#>    1.076399    1.652702 
inflation_factors(40, "S2")                       # an even n: interpolated quartiles, type 7
#>    kappa_mu kappa_sigma 
#>    1.070136    1.620424 
inflation_factors(21, "S4")                       # quartiles with the 10th and 90th percentiles
#>    kappa_mu kappa_sigma 
#>    1.028619    1.207533 
```
