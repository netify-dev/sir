# Confidence Intervals from Bootstrap SIR Results

Computes confidence intervals from a `boot_sir` object. For `block` and
`parametric` bootstraps these are percentile intervals from the
replicate distribution. For the `dyad` jackknife they are normal
(`estimate +/- z * se`) intervals from the jackknife standard error; the
jackknife produces a covariance, not a replicate distribution, so
percentiles do not apply.

## Usage

``` r
# S3 method for class 'boot_sir'
confint(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  A `boot_sir` object from
  [`boot_sir`](https://netify-dev.github.io/sir/reference/boot_sir.md).

- parm:

  Character vector of parameter names or integer indices. If NULL
  (default), returns intervals for all parameters.

- level:

  Confidence level between 0 and 1. Default is 0.95.

- ...:

  Additional arguments (unused).

## Value

A matrix with one row per parameter and columns for the lower and upper
bounds, labeled by percentage.
