# fivenumMeta 0.1.0

* First release: the functions that accompany the paper of Wu, Shi and Tong
  (2026), "Precision of the imputed mean and standard deviation from a
  five-number summary".

* `se_imputed_df()` takes a meta-analysis data set, one row per study, reads the
  reporting scenario off the columns each row has, and appends the imputed mean
  and SD with their naive and exact plug-in standard errors. A row with a sample
  mean and SD is the scenario `"sample"`, and a row with a mean, the minimum and
  the maximum and no SD is `"mean_range"`: the mean is kept, and the SD is imputed
  from the range. For two groups per row, `se_imputed_diff()` does the same and
  adds the effect of the study with its naive and its exact standard error, in the
  names `metafor::rma()` takes as they are: the difference of means
  (`measure = "MD"`, the default) or the log variability ratio (`measure = "lnVR"`).

* `pool_imputed()` pools the studies of `se_imputed_diff()` in one call, under the
  naive and under the exact standard errors side by side, by the common-effect
  model (`method = "EE"`, the default) or the DerSimonian-Laird random-effects
  model (`"DL"`): the pooled estimate, its standard error, confidence interval and
  p-value, Q, I^2, tau^2 and the weights. REML and the other estimators stay with
  metafor, which takes the same columns.

* `se_imputed()` returns one study's imputed mean and SD with their plug-in
  standard errors, as an object of class `"se_imputed"`: a list with a compact
  print. A row of `se_imputed_df()` whose standard errors would be zero (a
  sample SD of zero, a summary whose values are all equal) gets `NA` and is
  named in the warning with the other rows that cannot be used.

* `c_mu()`, `c_sigma()` and `inflation_factors()` give the exact factors at every
  sample size: a reported quantile is a fixed combination of two neighboring order
  statistics under the definition `type` (1 to 9 as in `quantile()`, default 7).
  So its moments, the estimators and their factors are exact at any n.
  `c_mu_shortcut()` and `c_sigma_shortcut()` give the closed-form shortcut formulas,
  for comparison and for use without the package.

* Every function defaults to the BLUE: `version = "blue"` in the study-level
  functions and `inflation_factors()`, `method = "blue"` in the others.
  `version = "wan_luo_shi"` and the `method` presets give the estimators of Wan,
  Luo and Shi, exact at every sample size as well.

* `estimator_weights()`, `scenario_ranks()` and `imputed_moments()` give the
  weights, the reported ranks and the second moments of the imputed statistics.

* The exact moments of the normal order statistics come from the package
  `osmoments`.

* The data set `vitd` carries the six studies of the vitamin D meta-analysis of
  Nnoaham and Clarke (2008), and the vignette `two-meta-analyses` reproduces the
  two meta-analyses of Wu, Shi and Tong (2026), one of means on `vitd` and one of
  SDs on `metamedian::dat.age`.
