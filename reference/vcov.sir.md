# Variance-Covariance Matrix for SIR Model Parameters

Returns the variance-covariance matrix of the estimated parameters. Two
types are available: classical (inverse Hessian) and robust (sandwich
estimator). The robust version is more reliable when the model is
misspecified or when observations are not independent.

## Usage

``` r
# S3 method for class 'sir'
vcov(object, type = c("classical", "robust"), ...)
```

## Arguments

- object:

  A fitted `sir` object from
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md).

- type:

  Character string: `"classical"` (default) for Hessian-based vcov, or
  `"robust"` for the sandwich estimator.

- ...:

  Additional arguments (unused).

## Value

A square matrix with rows and columns named by parameter. Returns NULL
if standard errors were not computed (`calc_se = FALSE`).

## See also

[`confint.sir`](https://netify-dev.github.io/sir/reference/confint.sir.md)
for confidence intervals.
