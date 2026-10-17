# Package index

## A study’s imputed mean and SD with their standard errors

The study-level function for practitioners.

- [`se_imputed()`](https://zwulab.github.io/fivenumMeta/reference/se_imputed.md)
  [`print(`*`<se_imputed>`*`)`](https://zwulab.github.io/fivenumMeta/reference/se_imputed.md)
  : Imputed mean and SD of one study, with their plug-in standard errors
- [`se_imputed_df()`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_df.md)
  : Imputed means and SDs with their standard errors for every study of
  a data frame
- [`se_imputed_diff()`](https://zwulab.github.io/fivenumMeta/reference/se_imputed_diff.md)
  : The effect of two groups per study, with its naive and exact
  standard errors

## Pooling

The pooled estimate under the naive and the exact standard errors, side
by side.

- [`pool_imputed()`](https://zwulab.github.io/fivenumMeta/reference/pool_imputed.md)
  [`print(`*`<pool_imputed>`*`)`](https://zwulab.github.io/fivenumMeta/reference/pool_imputed.md)
  : Pool the studies under the naive and the exact standard errors, side
  by side

## The exact factors

The standard deviation factors and the inflation factors of the imputed
mean and SD, under normality.

- [`c_mu()`](https://zwulab.github.io/fivenumMeta/reference/c_mu.md) :
  The scalar C_mu(n): standard deviation of the imputed mean in units of
  sigma
- [`c_sigma()`](https://zwulab.github.io/fivenumMeta/reference/c_sigma.md)
  : The scalar C_sigma(n): standard deviation of the imputed SD in units
  of sigma
- [`inflation_factors()`](https://zwulab.github.io/fivenumMeta/reference/inflation_factors.md)
  : Inflation factors of the imputed mean and SD
- [`imputed_moments()`](https://zwulab.github.io/fivenumMeta/reference/imputed_moments.md)
  : Second moments of the imputed mean and SD from a five-number summary

## The shortcut formulas

Closed forms in the sample size, with the constants as the companion
paper prints them.

- [`c_mu_shortcut()`](https://zwulab.github.io/fivenumMeta/reference/c_mu_shortcut.md)
  : Shortcut formula for C_mu(n)
- [`c_sigma_shortcut()`](https://zwulab.github.io/fivenumMeta/reference/c_sigma_shortcut.md)
  : Shortcut formula for C_sigma(n)

## Weights and ranks

- [`estimator_weights()`](https://zwulab.github.io/fivenumMeta/reference/estimator_weights.md)
  : Weight vectors of the published quantile-based estimators and the
  BLUE
- [`scenario_ranks()`](https://zwulab.github.io/fivenumMeta/reference/scenario_ranks.md)
  : Ranks of the order statistics reported under each scenario

## Data

The six studies of the vitamin D meta-analysis, the mean case of the
second vignette.

- [`vitd`](https://zwulab.github.io/fivenumMeta/reference/vitd.md) :
  Serum vitamin D in tuberculosis cases and controls: six studies

## Package

- [`fivenumMeta`](https://zwulab.github.io/fivenumMeta/reference/fivenumMeta-package.md)
  [`fivenumMeta-package`](https://zwulab.github.io/fivenumMeta/reference/fivenumMeta-package.md)
  : fivenumMeta: Five-Number Summaries in Meta-Analysis, Imputed Means
  and SDs with Standard Errors
