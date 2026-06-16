# Rolling-origin cross-validation for a fitted sir model

Performs expanding-window (rolling-origin) time-series cross-validation.
For each origin `o`, the model is refit on periods `1..o` and used to
forecast periods `o+1..o+horizon`, which are scored against the held-out
actuals. Training periods are always strictly before the test periods,
but forecasts are conditional on any future `Z` or dynamic-`W` slices
stored on `object`; interpret the result as conditional validation
unless those covariates are known, pre-specified, or generated as a
scenario at the forecast origin. A last-value-carried-forward naive
baseline is scored alongside the model so the result is a forecasting
horse race.

## Usage

``` r
cv_sir(
  object,
  initial = NULL,
  horizon = 1L,
  origins = NULL,
  baseline = TRUE,
  ...
)
```

## Arguments

- object:

  a fitted `sir`/`sir_fit` object; its data, family, and structural
  settings are reused for the refits.

- initial:

  integer; length of the first training window. Default
  `ceiling(n_periods / 2)`.

- horizon:

  integer; forecast horizon scored at each origin (default 1, one step
  ahead).

- origins:

  optional integer vector of training-window end points; defaults to
  `seq(initial, n_periods - horizon)`. Pass a sparse vector (e.g.
  `seq(20, 90, by = 5)`) to subsample origins and cut cost.

- baseline:

  logical; also score the last-value-carried-forward baseline (default
  `TRUE`).

- ...:

  passed to the
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md) refit. The
  structural settings reused from `object` — `family`, `method`,
  `fix_receiver`, `symmetric`, `bipartite`, `calc_se`, `W_recv` — are
  reserved: passing any of them here is ignored with a warning (every
  fold must fit the same model as `object`).

## Value

a `sir_cv` object with the per-origin score table, the aggregate (mean)
scores, the per-metric effective n (`eff_n`), the requested vs failed
fold counts, and (if requested) the baseline aggregate.

## Details

Each refit forces `calc_se = FALSE` (standard errors are not needed for
point forecasts and this is a large speed win). The cost is one model
refit per origin, and each refit's training window grows with the
origin, so the default (all origins from `initial` to
`n_periods - horizon`) is roughly quadratic in `n_periods` — on the
order of minutes for a hundred-period series. For long series subsample
`origins` (e.g. `seq(initial, n_periods - horizon, by = 5)`). The model
and the naive baseline are scored on the same observed cells each origin
so the comparison is fair; non-finite model predictions on those cells
make the fold fail rather than silently changing the denominator. With
`horizon > 1` the future cells are pooled across the horizon into one
score per origin. Folds that fail to fit, forecast, or score are dropped
with a warning and the aggregate averages over those that succeeded. For
a full-bilinear (`W_recv`) fit the refits draw random restarts; the
fit's `seed` is threaded through so the result is reproducible.

## See also

[`forecast.sir_fit`](https://netify-dev.github.io/sir/reference/forecast.sir_fit.md)
for the one-off forecasts CV scores,
[`score_sir`](https://netify-dev.github.io/sir/reference/score_sir.md)
for the scoring rules used.

## Examples

``` r
dat <- sim_sir(m = 10, T_len = 20, p = 2, q = 1, family = "poisson", seed = 1)
fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 1)
cv <- cv_sir(fit, initial = 12, origins = c(12, 15, 18))
cv
#> 
#> ── Rolling-origin cross-validation (3 origins, horizon 1) ──
#> 
#> Family: "poisson"
#> 
#> Model (out-of-sample, averaged over origins):
#> rmse: 1.101
#> mae: 0.83
#> deviance: 1.095
#> 
#> Naive (last value carried forward):
#> rmse: 1.693
#> mae: 1.185
#> deviance: 33.99
```
