# Weight vectors of the published quantile-based estimators and the BLUE

Returns the weight vector \\w\\ such that the estimator is \\w'Y\\, with
\\Y\\ the reported values in increasing order (for example the minimum,
the median and the maximum under S1). All formulas were read from the
source PDFs.

## Usage

``` r
estimator_weights(
  n,
  scenario = c("S3", "S1", "S2", "S4", "S5"),
  estimand = c("mean", "sd"),
  method = c("blue", "luo", "wan", "shi"),
  moments = NULL,
  type = 7
)
```

## Arguments

- n:

  A single integer sample size. The Wan, Luo and Shi presets are closed
  forms in \\n\\ and accept any \\n \ge 2\\; the BLUE needs the sample
  sizes of
  [`imputed_moments`](https://zwulab.github.io/fivenumMeta/reference/imputed_moments.md).

- scenario:

  Character. One of `"S1"` (min, median, max), `"S2"` (three quartiles),
  `"S3"` (all five). McGrath 2023 labels; note Wan 2014's C2/C3 are
  S3/S2 respectively. Also `"S4"` (the quartiles with the 10th and 90th
  percentiles) and `"S5"` (all five with the two percentiles), the names
  of Wu, Shi and Tong (2026) for Scenarios 5 and 7 of Balakrishnan et
  al. (2022); only the BLUE (`"blue"`) is defined for these two.

- estimand:

  `"mean"` or `"sd"`.

- method:

  One of the following; `"blue"` is the default. Under S4 and S5 only
  `"blue"`.

  `"blue"`

  :   The default. Lloyd's generalized-least-squares best linear
      unbiased estimator computed from the exact means and covariance
      matrix of the reported values (Balakrishnan 2022 eqs. 8-9; Yang
      2022 Definition 3.1).

  `"wan"`

  :   The estimators used/recommended in Wan et al. (2014). Mean: Hozo's
      \\(a+2m+b)/4\\ for S1 (Wan eq. 3), \\(q_1+m+q_3)/3\\ for S2 (Wan
      eq. 14), Bland's \\(a+2q_1+2m+2q_3+b)/8\\ for S3 (Wan eq. 10). SD:
      \\(b-a)/\xi_B(n)\\ for S1 (eq. 9), \\(q_3-q_1)/\eta_B(n)\\ for S2
      (eq. 16), the average of the two for S3 (eq. 13), with
      Blom-approximated divisors \\\xi_B(n) =
      2\Phi^{-1}\\(n-0.375)/(n+0.25)\\\\, \\\eta_B(n) =
      2\Phi^{-1}\\(0.75n-0.125)/(n+0.25)\\\\.

  `"luo"`

  :   Mean only. Luo et al. (2018) eqs. (7), (11), (15): S1 \\w_1 =
      4/(4+n^{0.75})\\ on the mid-range; S2 \\w = 0.7+0.39/n\\ on the
      mid-quartile range; S3 \\w_1 = 2.2/(2.2+n^{0.75})\\, \\w_2 =
      0.7-0.72/n^{0.55}\\.

  `"shi"`

  :   SD only, S3 only. Shi et al. (2020) eqs. (5), (9), (10): weight
      \\w = 1/(1+0.07 n^{0.6})\\ on \\(b-a)/\xi_B\\, rest on
      \\(q_3-q_1)/\eta_B\\. For S1/S2 Shi recommends Wan's estimators.

- moments:

  Optional precomputed result of `osmoments::os_moments(n, ranks)` for
  the standard normal, at this `n` and at the ranks that the reported
  values use (at `n = 4Q+1` under type 7 these are
  `scenario_ranks(n, scenario)`), or a list with its components `alpha`
  and `B` alone. The sample size, the distribution, the ranks, the
  location and scale and the shape of `alpha` and `B` are checked where
  the object carries them; a bare list is checked for its size only. The
  package keeps the moments it computes itself for the session, so
  repeated calls at one sample size integrate once without this
  argument.

- type:

  The definition of the sample quantiles, an integer from 1 to 9 as in
  [`quantile`](https://rdrr.io/r/stats/quantile.html); 7, R's default,
  unless a study states another.

## Value

Named numeric weight vector over the reported values, in increasing
order, named as the scenario reports them: `a` the minimum, `p10` the
10th percentile, `q1`, `m` the median, `q3`, `p90` the 90th percentile
and `b` the maximum.

## Details

The Wan/Luo/Shi weights are closed-form functions of \\n\\ and are
returned for any single integer \\n \ge 2\\, as the source papers apply
them. The `"blue"` weights are computed from the exact means and
covariance matrix of the reported values at any \\n\\ (see
[`imputed_moments`](https://zwulab.github.io/fivenumMeta/reference/imputed_moments.md),
argument `type`) and use `moments`.

## References

Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard
deviation from a five-number summary. Manuscript.

## Examples

``` r
estimator_weights(25, "S3", "mean")           # the BLUE of the mean, the default
#>          a         q1          m         q3          b 
#> 0.08214526 0.28565964 0.26439020 0.28565964 0.08214526 
estimator_weights(25, "S3", "sd")             # the BLUE of the SD
#>             a            q1             m            q3             b 
#> -1.712665e-01 -2.565653e-01 -2.315005e-16  2.565653e-01  1.712665e-01 
estimator_weights(25, "S3", "mean", "luo")
#>          a         q1          m         q3          b 
#> 0.08221017 0.28870353 0.25817261 0.28870353 0.08221017 
estimator_weights(25, "S3", "sd", "shi")
#>          a         q1          m         q3          b 
#> -0.1716592 -0.2559979  0.0000000  0.2559979  0.1716592 
```
