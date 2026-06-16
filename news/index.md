# Changelog

## sir 1.0.0

- `sir(symmetric = TRUE)` fits a genuine undirected model with a single
  shared influence operator `A = B = sum_k gamma_k W_k` (the quadratic
  form `A X A'`). All `gamma_k` are estimated; the operator is
  identified up to global sign, fixed so the largest-magnitude `gamma_k`
  is positive. `sim_sir(symmetric = TRUE)` simulates a matching process.
  (`symmetric` and `fix_receiver` are mutually exclusive.)

- Standard errors default to an actor-clustered cluster-robust sandwich
  for [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html), and
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) (HC1
  factor, `t(G - 1)` reference), for both directed and symmetric fits.
  [`summary()`](https://rdrr.io/r/base/summary.html) still prints the
  classical SE; request it anywhere with `se.type = "classical"`. For
  dynamic (4D) `W` the default accessors fall back to classical Wald;
  use `boot_sir(type = "dyad")` there.

- [`boot_sir()`](https://netify-dev.github.io/sir/reference/boot_sir.md)
  supports symmetric fits for all bootstrap types, bipartite cluster SEs
  count senders and receivers as distinct actors, and
  [`predict()`](https://rdrr.io/r/stats/predict.html) errors on
  covariate arrays passed outside `newdata`.

- [`sim_sir()`](https://netify-dev.github.io/sir/reference/sim_sir.md)
  warns when user-supplied influence coefficients imply a non-stationary
  (spectral gain \>= 1) Poisson/Normal process, instead of silently
  returning an explosive series.
