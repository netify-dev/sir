# Bayesian Information Criterion for a SIR Model

Computes BIC = -2 \* log-likelihood + log(nobs) \* (number of
parameters). BIC penalizes model complexity more heavily than AIC for
large samples.

## Usage

``` r
# S3 method for class 'sir'
BIC(object, ...)
```

## Arguments

- object:

  A fitted `sir` object from
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md).

- ...:

  Additional arguments for comparison with other models.

## Value

Numeric BIC value. Lower is better.
