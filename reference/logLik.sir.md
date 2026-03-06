# Extract Log-Likelihood from a SIR Model

Returns the log-likelihood at convergence as a `logLik` object, with
attributes for degrees of freedom and number of observations. This
allows [`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) to work directly on the
result.

## Usage

``` r
# S3 method for class 'sir'
logLik(object, ...)
```

## Arguments

- object:

  A fitted `sir` object from
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md).

- ...:

  Additional arguments (unused).

## Value

A `logLik` object with attributes `df` (number of estimated parameters)
and `nobs` (number of observations).
