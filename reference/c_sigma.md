# The scalar C_sigma(n): standard deviation of the imputed SD in units of sigma

\\C\_\sigma(n) = \sqrt{w\_\sigma' B w\_\sigma}\\ so that
\\\mathrm{SD}(\hat\sigma) = C\_\sigma(n)\\\sigma\\, the true standard
deviation of the imputed SD (a function of the unknown population
\\\sigma\\); the plug-in standard error is \\\mathrm{SE}(\hat\sigma) =
C\_\sigma(n)\\\hat\sigma\\. The factor is exact (numerical integration
of the order-statistic moments), and
[`c_sigma_shortcut`](https://zwulab.github.io/fivenumMeta/reference/c_sigma_shortcut.md)
gives its approximation formula. This is the SD companion of
[`c_mu`](https://zwulab.github.io/fivenumMeta/reference/c_mu.md).

## Usage

``` r
c_sigma(
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

  SD estimator preset (see
  [`estimator_weights`](https://zwulab.github.io/fivenumMeta/reference/estimator_weights.md)).
  The default `"blue"` is the BLUE of \\\sigma\\, the default estimator
  version of Wu, Shi and Tong (2026); `NULL` selects the SD estimator of
  Wan (S1 and S2) or Shi (S3).

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

Numeric scalar \\C\_\sigma(n)\\.

## Details

Assumptions as in
[`imputed_moments`](https://zwulab.github.io/fivenumMeta/reference/imputed_moments.md):
iid normal data, the reported values following definition `type`. The
value is `sqrt(imputed_moments(...)$var_sd)`.

## References

Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard
deviation from a five-number summary. Manuscript.

## Examples

``` r
c_sigma(25, "S2")                    # exact C_sigma(n) for the BLUE of sigma under S2
#> [1] 0.237215
c_sigma(25, "S2", method = NULL)     # the same for Wan's estimator
#> [1] 0.2375384
c_sigma(25, "S2") * sqrt(2 * 25)     # inflation over the large-sample sigma/sqrt(2n)
#> [1] 1.677364
```
