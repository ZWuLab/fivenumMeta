# Shortcut formula for C_sigma(n)

Closed-form shortcut formula for the exact
[`c_sigma`](https://zwulab.github.io/fivenumMeta/reference/c_sigma.md).
The form differs by scenario for a structural reason. Under S2 and S3
\\n\\C^2\_\sigma(n)\\ has the finite limit \\h\\ (Theorem 1 of Wu, Shi
and Tong, 2026), so the formula is \\A/\sqrt{n + \dots}\\ with \\A =
\sqrt{h}\\ rounded to four decimals. Under S3 the minimum and maximum
add precision of order \\(\log n)^2\\, against order \\n\\ from the
quartiles, which the terms \\c_1(\log n)^2 + c_2\log n\\ stand for.
Under S1 the range-based factor has no such limit and decays only like
\\1/\log n\\, so the formula is the logarithmic \\c_1/(\log n +
c_2)^{c_3}\\. Two constant sets are carried: `version = "blue"` (the
default) for the factor of the BLUE of \\\sigma\\, and
`version = "wan_luo_shi"` for the factor of the SD estimators of Wan
(under S1 and S2) and Shi (under S3). Both are fitted by minimax
relative error over the paper's 39-point grid of \\n = 4Q+1\\ (every
such \\n\\ from 5 to 101 and 14 larger ones up to 1001), with the
constants exactly as the paper prints them.

## Usage

``` r
c_sigma_shortcut(
  n,
  scenario = c("S3", "S1", "S2"),
  version = c("blue", "wan_luo_shi")
)
```

## Arguments

- n:

  Integer sample size(s); vectorized. Applied at any \\n\\, as the
  source papers apply their formulas. No upper bound on \\n\\ is
  checked; see Details for the accuracy beyond the grid.

- scenario:

  One of "S1", "S2", "S3".

- version:

  `"blue"` (default) or `"wan_luo_shi"`.

## Value

Numeric vector of approximate \\C\_\sigma(n)\\ values.

## Details

Accuracy is the relative error of the printed constants against the
exact factor, on the grid and at \\n = 2001, 5001, 10001\\ beyond it.
BLUE:

- S1: \\0.8278/(\log n + 0.330)^{1.2041}\\; error \\\le 0.37\\\\ on the
  grid, and \\-0.8\\\\, \\-1.4\\\\, \\-1.9\\\\ at \\n = 2001, 5001,
  10001\\: this formula has no limit to hold it, and its error grows
  slowly beyond the grid.

- S2: \\1.1664/\sqrt{n - 0.863}\\; error \\\le 0.13\\\\ on the grid and
  \\\le 0.003\\\\ beyond.

- S3: \\1.1664/\sqrt{n + 5.267(\log n)^2 - 10.364\log n + 8.327}\\;
  error \\\le 0.42\\\\ on the grid and \\\le 0.32\\\\ beyond.

Shortcut version:

- S1: \\0.8705/(\log n + 0.418)^{1.2199}\\; error \\\le 0.37\\\\ on the
  grid, and \\-0.8\\\\, \\-1.5\\\\, \\-2.0\\\\ beyond.

- S2: \\1.1664/\sqrt{n - 0.845}\\; error \\\le 0.099\\\\ on the grid and
  \\\le 0.002\\\\ beyond.

- S3: \\1.1664/\sqrt{n + 5.191(\log n)^2 - 10.204\log n + 8.525}\\;
  error \\\le 0.44\\\\ on the grid and \\\le 0.27\\\\ beyond.

Assumptions: iid normal data. The errors are measured on the grid of \\n
= 4Q+1\\, where the quartiles are single order statistics (type 7), and
at the three larger sizes named; they are not bounds for every \\n\\.
[`c_sigma`](https://zwulab.github.io/fivenumMeta/reference/c_sigma.md)
computes the exact factor at any \\n\\ and definition.

## References

Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard
deviation from a five-number summary. Manuscript.

## Examples

``` r
c_sigma_shortcut(25, "S2")
#> [1] 0.2374137
c_sigma_shortcut(5000, "S3")    # no upper limit on n
#> [1] 0.01601851
```
