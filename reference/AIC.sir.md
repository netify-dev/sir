# Akaike Information Criterion for a SIR Model

Computes AIC = -2 \* log-likelihood + k \* (number of parameters). Use
this to compare SIR models with different specifications (e.g.,
different numbers of influence covariates).

## Usage

``` r
# S3 method for class 'sir'
AIC(object, ..., k = 2)
```

## Arguments

- object:

  A fitted `sir` object from
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md).

- ...:

  Additional arguments for comparison with other models.

- k:

  Numeric penalty per parameter (default 2 for standard AIC).

## Value

Numeric AIC value. Lower is better.
