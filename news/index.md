# Changelog

## sir 1.1.0

- [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html), and
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) now report
  the actor-clustered cluster-robust sandwich for dynamic (4D) `W` fits
  too (directed, symmetric, and `fix_receiver`), with the same
  `t(G - 1)` reference as the static case. Previously these fell back to
  classical Wald intervals for dynamic `W`.

- Standard errors are always available now. When an analytic covariance
  cannot be formed (full-bilinear bipartite fits, or an ill-conditioned
  Hessian) and the model converged,
  [`sir()`](https://netify-dev.github.io/sir/reference/sir.md)
  automatically attaches a delete-one-actor jackknife covariance, so
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`summary()`](https://rdrr.io/r/base/summary.html), and
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) keep
  working without a manual
  [`boot_sir()`](https://netify-dev.github.io/sir/reference/boot_sir.md)
  call. The fit records `se_source = "jackknife"` when this fallback
  fires (`NULL` otherwise), and the printouts label which kind of
  standard error they are showing.

- Fixed a crash when fitting dynamic (4D) `W` with a single influence
  covariate (`p = 1`): the per-period slice `W[, , , t]` collapsed to a
  matrix and failed the C++ cube conversion in the operator and
  fitted-value construction. Single-covariate dynamic models now fit.

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
  classical SE; request it anywhere with `se.type = "classical"`.

- [`boot_sir()`](https://netify-dev.github.io/sir/reference/boot_sir.md)
  supports symmetric fits for all bootstrap types, bipartite cluster SEs
  count senders and receivers as distinct actors, and
  [`predict()`](https://rdrr.io/r/stats/predict.html) errors on
  covariate arrays passed outside `newdata`.

- [`sim_sir()`](https://netify-dev.github.io/sir/reference/sim_sir.md)
  warns when user-supplied influence coefficients imply a non-stationary
  (spectral gain \>= 1) Poisson/Normal process, instead of silently
  returning an explosive series.
