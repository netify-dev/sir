# Predictions from a Fitted SIR Model

Generates predictions from a fitted SIR model for the training data or
for new data. Predictions can be on the link scale (linear predictor) or
the response scale (expected counts, probabilities, or means).

## Usage

``` r
# S3 method for class 'sir'
predict(object, newdata = NULL, type = c("response", "link"), ...)
```

## Arguments

- object:

  A fitted `sir` object from
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md).

- newdata:

  Optional named list with components `W` (3D or 4D array; 3D only for
  full-bilinear bipartite fits), `X` (3D array), `Z` (3D or 4D array),
  and for full-bilinear bipartite fits `W_recv` (3D receiver-side
  array), for scenario prediction. Dimensions must match the original
  fit. Any component not supplied is taken from the original fit. All
  supplied time-varying components must agree on the number of time
  periods, and `Z` must carry the same number of covariates as the fit;
  a genuine mismatch is an error rather than being silently recycled. If
  NULL (default), returns predictions for the training data. Note:
  unlike many R predict methods, `newdata` is a list of arrays, not a
  data frame. For a full-bilinear bipartite fit, supply `newdata$W_recv`
  to vary the receiver-side influence structure; otherwise the fitted
  `W_recv` is reused. Unlike the square one-mode case the bipartite
  diagonal is a genuine prediction and is not set to NA.

- type:

  Character string: `"link"` for linear predictor (eta) or `"response"`
  for expected values on the original scale. Default is `"response"`.

- ...:

  Additional arguments (unused).

## Value

An array (n1 x n2 x T) of predicted values on the requested scale.

## Details

For model-implied scenario analysis, supply modified arrays in
`newdata`. For example, to see how fitted values change when a covariate
is increased by one unit, pass the modified Z array while keeping W and
X from the original fit. Causal counterfactual interpretation requires
additional design assumptions.

## Examples

``` r
dat <- sim_sir(m = 10, T_len = 20, p = 2, q = 1, family = "poisson", seed = 1)
fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 1)
# in-sample fitted values (response scale)
pred <- predict(fit)
# scenario: only W/X/Z are read from newdata (Y is ignored)
Zcf <- dat$Z; Zcf[, , 1, ] <- Zcf[, , 1, ] + 1
pred_cf <- predict(fit, newdata = list(W = dat$W, X = dat$X, Z = Zcf))
```
