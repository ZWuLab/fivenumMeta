# Second moments of the imputed mean and SD from a five-number summary

Computes \\\mathrm{Var}(\hat\mu)/\sigma^2\\,
\\\mathrm{Var}(\hat\sigma)/\sigma^2\\,
\\\mathrm{Cov}(\hat\mu,\hat\sigma)/\sigma^2\\ and \\C(n) =
\mathrm{SD}(\hat\mu)/\sigma\\ for linear estimators \\\hat\mu = w\_\mu'
Y\\, \\\hat\sigma = w\_\sigma' Y\\ of the reported values \\Y\\, under
normality, using exact order-statistic moments. The identity used is
\\\mathrm{Var}(w'Y) = \sigma^2 w'B_Yw\\ (Balakrishnan et al. 2022, eqs.
49-51), where \\B_Y\\ is the covariance matrix of the standardized
reported values. The default, `method_mean = "blue"` with
`method_sd = "blue"`, is the BLUE, as in every function of the package;
`method_mean = "luo"` with `method_sd = NULL` gives the estimators of
Wan, Luo and Shi.

## Usage

``` r
imputed_moments(
  n,
  scenario = c("S3", "S1", "S2", "S4", "S5"),
  method_mean = "blue",
  method_sd = "blue",
  w_mean = NULL,
  w_sd = NULL,
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

- method_mean:

  Character, NULL or NA. One of `"blue"` (default), `"luo"` or `"wan"`,
  the mean estimator presets of
  [`estimator_weights`](https://zwulab.github.io/fivenumMeta/reference/estimator_weights.md).
  Ignored if `w_mean` is given. `NULL` or `NA` skips the mean part.

- method_sd:

  Character, NULL or NA. One of `"blue"` (default), `"wan"`, `"shi"` (S3
  only). `NULL` pairs Wan's SD under S1 and S2 with Shi's under S3, the
  pairing of the Wan-Luo-Shi version (Shi 2020 sec. 2 recommends Wan's
  SD for S1/S2). Ignored if `w_sd` given. Set to `NA` to skip the SD
  part.

- w_mean, w_sd:

  Optional numeric weight vectors overriding the presets (length =
  number of reported values, in increasing order).

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

A list with components: `var_mean` \\= w\_\mu'Bw\_\mu\\, `var_sd` \\=
w\_\sigma'Bw\_\sigma\\, `cov` \\= w\_\mu'Bw\_\sigma\\ (all in units of
\\\sigma^2\\); `c_n` \\= \sqrt{w\_\mu'Bw\_\mu}\\; `w_mean`, `w_sd`;
`gaps`: named vector of unbiasedness gaps (`mean_sum` \\= w\_\mu'1 -
1\\, `mean_alpha` \\= w\_\mu'\alpha\\ = bias of \\\hat\mu\\ in
\\\sigma\\ units, `sd_sum` \\= w\_\sigma'1\\, `sd_alpha` \\=
w\_\sigma'\alpha - 1\\ = relative bias of \\\hat\sigma\\); `ranks` (the
order statistics the reported values use), `n`, `scenario`, `type`.

## Details

All second moments are proportional to the true \\\sigma^2\\; nothing
here plugs in an estimate of \\\sigma\\. In particular
\\\mathrm{SD}(\hat\mu) = C(n)\sigma\\, and the plug-in standard error is
\\C(n)\hat\sigma\\, not \\\hat\sigma/\sqrt{n}\\. Each reported value is
an order statistic (the minimum, the maximum) or a sample quantile,
which under every definition of
[`quantile`](https://rdrr.io/r/stats/quantile.html) is a fixed linear
combination \\(1-h)X\_{(j)} + hX\_{(j+1)}\\ of two adjacent order
statistics (Hyndman and Fan 1996). So \\Y = W'X\\, with \\W\\ from
[`osmoments::os_quantile_weights()`](https://zwulab.github.io/osmoments/reference/os_quantile_weights.html),
and the means and covariance matrix of the standardized \\Y\\ are
\\\alpha_Y = W'\alpha\\ and \\B_Y = W'BW\\, exact at every \\n\\. At \\n
= 4Q+1\\ under type 7 the quartiles are the order statistics at ranks
\\Q+1, 2Q+1, 3Q+1\\; at other \\n\\, or under another definition, they
are interpolated, and the result depends on the definition. Assumptions:
i.i.d. normal sample; the reported values follow definition `type`.

## References

Balakrishnan N, Rychtar J, Taylor D, Walter SD (2022). Unified approach
to optimal estimation of mean and standard deviation from sample
summaries. Statistical Methods in Medical Research 31(11), 2087-2103.
(eqs. 49-51) Luo D, Wan X, Liu J, Tong T (2018). Optimally estimating
the sample mean from the sample size, median, mid-range, and/or
mid-quartile range. Statistical Methods in Medical Research 27(6),
1785-1805. Shi J, Luo D, Weng H, Zeng XT, Lin L, Chu H, Tong T (2020).
Optimally estimating the sample standard deviation from the five-number
summary. Research Synthesis Methods 11(5), 641-654. Wan X, Wang W, Liu
J, Tong T (2014). Estimating the sample mean and standard deviation from
the sample size, median, range and/or interquartile range. BMC Medical
Research Methodology 14, 135. Yang X, Hutson AD, Wang D (2022). A
generalized BLUE approach for combining location and scale information
in a meta-analysis. Journal of Applied Statistics 49(15), 3846-3867.
(eq. 22) Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and
standard deviation from a five-number summary. Manuscript.

## Examples

``` r
im <- imputed_moments(25, "S3")            # the BLUE of the mean and of the SD, the default
im$c_n                                     # SD(mu-hat)/sigma
#> [1] 0.2068675
imputed_moments(25, "S3", "luo", NULL)$c_n # Luo's mean, with Shi's SD
#> [1] 0.2068703
sqrt(1/25)                                 # naive sigma/sqrt(n) factor, for contrast
#> [1] 0.2
```
