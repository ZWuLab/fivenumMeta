# Shortcut formula for C_mu(n)

Closed-form shortcut formula \\C\_\mu(n) \approx A/\sqrt{n + c_1 \log
n + c_2}\\ for the exact
[`c_mu`](https://zwulab.github.io/fivenumMeta/reference/c_mu.md), in the
style of Wan's \\\xi_B(n)\\ and Luo's \\w\_{opt}(n)\\ formulas. The
leading constant \\A\\ is the limit of \\\sqrt{n}\\C\_\mu(n)\\ (Theorem
1 of Wu, Shi and Tong, 2026) rounded to four decimals: \\\sqrt{\pi/2}\\
under S1 and \\\sqrt{v}\\ under S2 and S3. The term \\c_1 \log n\\
stands for the precision that the minimum and maximum add, which is of
order \\\log n\\ against order \\n\\ from the central order statistics;
under S2 no extreme is reported and \\c_1 = 0\\. This reason is the
BLUE's; for the Wan-Luo-Shi version the same form is an empirical fit.
Two constant sets are carried: `version = "blue"` (the default) for the
factor of the BLUE, the default estimator version of Wu, Shi and Tong
(2026), and `version = "wan_luo_shi"` for the factor of Luo's mean
estimator. Both are fitted by minimax relative error over the paper's
39-point grid of \\n = 4Q+1\\ (every such \\n\\ from 5 to 101 and 14
larger ones up to 1001), and are used here exactly as the paper prints
them (A to four decimals, the c's to three), so that the software and
the paper carry the same formula. The errors below are those of the
printed constants, on the grid and at \\n = 2001, 5001, 10001\\ beyond
it, computed by the fitting script. The formulas serve a reader without
this package and a comparison with the exact values;
[`c_mu`](https://zwulab.github.io/fivenumMeta/reference/c_mu.md)
computes the exact factor at every \\n\\, and
[`se_imputed`](https://zwulab.github.io/fivenumMeta/reference/se_imputed.md)
uses it.

## Usage

``` r
c_mu_shortcut(
  n,
  scenario = c("S3", "S1", "S2"),
  version = c("blue", "wan_luo_shi")
)
```

## Arguments

- n:

  Integer sample size(s); vectorized. The formulas are continuous in
  \\n\\ and are applied at any \\n\\, as the source papers apply theirs.
  Each tends to its limit as \\n\\ grows, so no upper bound on \\n\\ is
  checked.

- scenario:

  One of "S1", "S2", "S3".

- version:

  `"blue"` (default) or `"wan_luo_shi"`: which estimator version's
  factor the formula approximates.

## Value

Numeric vector of approximate \\C\_\mu(n)\\ values.

## Details

Fitted forms and accuracy (relative error against the exact factor),
BLUE:

- S1: \\1.2533/\sqrt{n + 2.006\log n - 0.810}\\; error \\\le 0.34\\\\ on
  the grid and \\\le 0.08\\\\ at \\n = 2001, 5001, 10001\\.

- S2: \\1.0780/\sqrt{n + 0.142}\\; error \\\le 0.046\\\\ on the grid and
  \\\le 0.003\\\\ beyond.

- S3: \\1.0780/\sqrt{n + 0.838\log n - 0.576}\\; error \\\le 0.33\\\\ on
  the grid and \\\le 0.09\\\\ beyond.

Shortcut version (Luo's shortcut mean):

- S1: \\1.2533/\sqrt{n + 2.006\log n - 0.811}\\; error \\\le 0.34\\\\ on
  the grid and \\\le 0.06\\\\ beyond.

- S2: \\1.0780/\sqrt{n + 0.142}\\; error \\\le 0.046\\\\ on the grid and
  \\\le 0.003\\\\ beyond.

- S3: \\1.0780/\sqrt{n + 0.835\log n - 0.572}\\; error \\\le 0.34\\\\ on
  the grid and \\\le 0.08\\\\ beyond.

Assumptions: iid normal data. The errors are measured on the grid of \\n
= 4Q+1\\, where the quartiles are single order statistics (type 7), and
at the three larger sizes named; they are not bounds for every \\n\\.
[`c_mu`](https://zwulab.github.io/fivenumMeta/reference/c_mu.md)
computes the exact factor at any \\n\\ and definition.

## References

Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard
deviation from a five-number summary. Manuscript.

## Examples

``` r
c_mu_shortcut(25, "S2")
#> [1] 0.2149903
c_mu_shortcut(25, "S2", version = "wan_luo_shi")
#> [1] 0.2149903
c_mu_shortcut(5000, "S3")        # no upper limit on n
#> [1] 0.01523523
```
