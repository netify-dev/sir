# Extract Model Coefficients from a SIR Model

Returns the estimated parameter values from a fitted SIR model. These
include exogenous covariate effects (theta), sender influence weights
(alpha), and receiver influence weights (beta).

## Usage

``` r
# S3 method for class 'sir'
coef(object, ...)
```

## Arguments

- object:

  A fitted `sir` object from
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md).

- ...:

  Additional arguments (unused).

## Value

Named numeric vector of estimated coefficients, tagged as in
[`summary()`](https://rdrr.io/r/base/summary.html): `(Z)` = exogenous
direct effects (theta); `(alphaW)` = sender-influence weights
(alpha_2..alpha_p for a directed fit, since alpha_1 = 1 is fixed);
`(betaW)` = receiver-influence weights (beta); `(gammaW)` = shared
undirected weights for a symmetric fit (all gamma_1..gamma_p estimated,
identified up to sign).

## See also

[`confint.sir`](https://netify-dev.github.io/sir/reference/confint.sir.md)
for confidence intervals,
[`vcov.sir`](https://netify-dev.github.io/sir/reference/vcov.sir.md) for
the variance-covariance matrix.
