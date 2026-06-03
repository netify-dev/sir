# Forecast future networks from a fitted sir model

Produces out-of-sample expected networks `h` steps beyond the end of the
estimation sample. Forecasts are iterated (plug-in): each period's
predicted mean is fed forward as the next-period lagged signal `X`,
transformed the same way
[`sim_sir`](https://netify-dev.github.io/sir/reference/sim_sir.md)
builds `X` from a lagged outcome (`log(Y + 1) / infl_scale` for
`poisson`, raw `Y / infl_scale` otherwise). This assumes `X` was built
with the same `infl_scale` used below; if it was not, supply
`infl_scale` — the function warns when the stored `X` disagrees.
Forecasts are conditional on any supplied future covariates, so
`Z_future` and `W_future` should be values that are known,
pre-specified, or generated as a scenario at the forecast origin.

## Usage

``` r
# S3 method for class 'sir_fit'
forecast(
  object,
  h = 1L,
  Z_future = NULL,
  W_future = NULL,
  Y_last = NULL,
  infl_scale = NULL,
  ...
)

# S3 method for class 'sir'
forecast(
  object,
  h = 1L,
  Z_future = NULL,
  W_future = NULL,
  Y_last = NULL,
  infl_scale = NULL,
  ...
)
```

## Arguments

- object:

  a fitted `sir`/`sir_fit` object.

- h:

  integer \>= 1; forecast horizon (number of future periods).

- Z_future:

  exogenous covariates for the future periods. Required when the model
  has `q > 0` predictors. These are conditioning values, not forecast by
  [`forecast()`](https://generics.r-lib.org/reference/forecast.html):
  pass only covariates known or fixed at the forecast origin, or values
  from an explicit scenario. Either an `n1 x n2 x q x h` array, or (when
  `h == 1`) an `n1 x n2 x q` array. Ignored with a message when
  `q == 0`.

- W_future:

  future influence covariates; required only when `object$dynamic_W` is
  `TRUE` and the model has influence covariates. As with `Z_future`,
  pass known, fixed, or scenario values. An `n1 x n1 x p x h` array (or
  `n1 x n1 x p` when `h == 1`). Static-W models reuse `object$W`.

- Y_last:

  optional `n1 x n2` matrix giving the most recent observed outcome that
  seeds the first forecast lag. Defaults to the last period of
  `object$Y`. Missing (NA) cells are treated as a zero lag.

- infl_scale:

  optional positive scalar; the divisor applied to the lagged outcome
  when building the forecast `X`. Defaults to `max(m - 1, 1)` for a
  one-mode network, `max(n1 - 1, 1)` for a sender-side bipartite fit,
  and `sqrt((n1 - 1)(n2 - 1))` for a full-bilinear bipartite fit.
  Override it when `X` was constructed with a different scaling.

- ...:

  unused.

## Value

an `n1 x n2 x h` array of expected outcomes on the response scale, with
actor dimnames carried from `object$Y` when present and a third
dimension labelled `h1, h2, ...`. Carries attribute `"family"`.

## Details

For `h == 1` the forecast equals
`predict(object, newdata = list(X = X_next, ...))` where `X_next` is the
transform of `Y_last` (off the diagonal for a square one-mode network,
whose self-ties are NA in the forecast but computed by `predict`). For
`h > 1` the period-`s` predicted mean becomes the lag for period
`s + 1`. This is a plug-in point forecast: it feeds the conditional mean
forward and does not propagate forecast uncertainty, so it produces no
predictive interval; for a nonlinear link (`poisson`/`binomial`) the
multi-step path is a biased approximation of the true conditional mean
(Jensen's inequality), growing with the horizon. The diagonal
(self-ties) is not a meaningful forecast for a square one-mode network
and is returned as `NA`.

## Examples

``` r
dat <- sim_sir(m = 10, T_len = 20, p = 2, q = 1, family = "poisson", seed = 1)
train_t <- 19
fit <- sir(dat$Y[, , 1:train_t], W = dat$W,
           X = dat$X[, , 1:train_t],
           Z = dat$Z[, , , 1:train_t],
           family = "poisson", fix_receiver = TRUE, seed = 1)
# one-step-ahead forecast, conditional on pre-specified future Z.
# In this simulated example we reuse the held-out Z slice to show the API.
Z_next <- dat$Z[, , , train_t + 1, drop = FALSE]
fc1 <- forecast(fit, h = 1, Z_future = Z_next)
dim(fc1)
#> [1] 10 10  1
score_sir(dat$Y[, , train_t + 1, drop = FALSE], fc1, "poisson")
#>      rmse       mae  deviance 
#> 1.0185078 0.7974237 1.0695265 
```
