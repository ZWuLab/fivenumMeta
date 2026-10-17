# The scalar C_mu(n): standard deviation of the imputed mean in units of sigma

\\C\_\mu(n) = \sqrt{w\_\mu' B w\_\mu}\\ so that \\\mathrm{SD}(\hat\mu) =
C\_\mu(n)\\\sigma\\ (the true standard deviation of the estimator, a
function of the unknown population \\\sigma\\). The plug-in standard
error is \\\mathrm{SE}(\hat\mu) = C\_\mu(n)\\\hat\sigma\\; its reference
distribution is the subject of a separate study. Exact (numerical
integration); for the closed-form shortcut formula see
[`c_mu_shortcut`](https://zwulab.github.io/fivenumMeta/reference/c_mu_shortcut.md).

## Usage

``` r
c_mu(
  n,
  scenario = c("S3", "S1", "S2", "S4", "S5"),
  method = "blue",
  moments = NULL,
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

- method:

  Mean estimator preset; see
  [`estimator_weights`](https://zwulab.github.io/fivenumMeta/reference/estimator_weights.md).
  The default `"blue"` is the BLUE, the default estimator version of Wu,
  Shi and Tong (2026); `"luo"` gives the factor of Luo's shortcut mean
  estimator.

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

Numeric scalar \\C(n)\\.

## References

Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard
deviation from a five-number summary. Manuscript.

## Examples

``` r
c_mu(25, "S2")                 # exact C(n) for the BLUE of the mean under S2
#> [1] 0.2150631
c_mu(25, "S2", "luo")          # the same for Luo's estimator
#> [1] 0.2150631
c_mu(25, "S2") * sqrt(25)      # inflation factor over sigma/sqrt(n)
#> [1] 1.075315
```
