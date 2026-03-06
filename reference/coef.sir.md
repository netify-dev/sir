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

Named numeric vector of estimated coefficients.

## See also

[`confint.sir`](https://netify-dev.github.io/sir/reference/confint.sir.md)
for confidence intervals,
[`vcov.sir`](https://netify-dev.github.io/sir/reference/vcov.sir.md) for
the variance-covariance matrix.
