# Ranks of the order statistics reported under each scenario

For \\n = 4Q+1\\: S1 reports ranks \\(1, 2Q+1, 4Q+1)\\ (min, median,
max), S2 ranks \\(Q+1, 2Q+1, 3Q+1)\\ (the quartiles), S3 all five (Luo
2018 sec. 2.2; McGrath 2023 scenario labels). For \\n = 20m+1\\: S4
reports the quartiles with the 10th and 90th percentiles at ranks
\\2m+1\\ and \\18m+1\\, that is \\p(n-1)+1\\ for \\p = 0.1, 0.9\\, and
S5 adds the minimum and maximum (Balakrishnan 2022 sec. 3.1, Scenarios 5
and 7).

## Usage

``` r
scenario_ranks(n, scenario = c("S3", "S1", "S2", "S4", "S5"))
```

## Arguments

- n:

  A single integer sample size at which every reported value is one
  order statistic under type 7: odd \\n \ge 3\\ under S1, \\n = 4Q + 1\\
  with \\Q \ge 1\\ under S2 and S3, and \\n = 20m + 1\\ with \\m \ge 1\\
  under S4 and S5. Any other \\n\\ is an error here, while the other
  functions of the package accept it.

- scenario:

  Character. One of `"S1"` (min, median, max), `"S2"` (three quartiles),
  `"S3"` (all five). McGrath 2023 labels; note Wan 2014's C2/C3 are
  S3/S2 respectively. Also `"S4"` (the quartiles with the 10th and 90th
  percentiles) and `"S5"` (all five with the two percentiles), the names
  of Wu, Shi and Tong (2026) for Scenarios 5 and 7 of Balakrishnan et
  al. (2022); only the BLUE (`"blue"`) is defined for these two.

## Value

Integer vector of ranks, increasing.

## Details

These are the sample sizes at which every reported value is a single
order statistic under R's default quantile definition (type 7). At other
sample sizes the functions of the package build the reported values from
two adjacent order statistics each (see
[`imputed_moments`](https://zwulab.github.io/fivenumMeta/reference/imputed_moments.md)).

## References

Wu Z, Shi J, Tong T (2026). Precision of the imputed mean and standard
deviation from a five-number summary. Manuscript.

## Examples

``` r
scenario_ranks(25, "S3")
#> [1]  1  7 13 19 25
```
