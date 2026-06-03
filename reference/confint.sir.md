# Confidence Intervals for SIR Model Parameters

Computes confidence intervals using either Wald-based intervals
(cluster- robust by default for supported static fits) or bootstrap
percentile intervals. Wald intervals require that the model was fit with
`calc_se = TRUE`. Bootstrap intervals require a
[`boot_sir`](https://netify-dev.github.io/sir/reference/boot_sir.md)
result and tend to be more reliable when the Hessian is ill-conditioned.

## Usage

``` r
# S3 method for class 'sir'
confint(
  object,
  parm = NULL,
  level = 0.95,
  boot = NULL,
  se.type = c("cluster", "classical", "robust", "twoway", "dyad"),
  ...
)
```

## Arguments

- object:

  A fitted `sir` object from
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md).

- parm:

  Character vector of parameter names or integer indices to include. If
  NULL (default), returns intervals for all parameters.

- level:

  Confidence level between 0 and 1. Default is 0.95.

- boot:

  Optional `boot_sir` object from
  [`boot_sir`](https://netify-dev.github.io/sir/reference/boot_sir.md).
  When provided, percentile intervals from the bootstrap distribution
  are used instead of Wald intervals.

- se.type:

  Character. For Wald intervals, which standard errors to use:
  `"cluster"` (default; multiway sender, receiver, and time clustering
  for supported static fits), `"classical"` (inverse-Hessian),
  `"robust"` (HC0 sandwich), `"twoway"` (alias for `"cluster"`), or
  `"dyad"` (directed-dyad clustering). Cluster intervals still rely on
  the Hessian bread being stable. Ignored when `boot` is supplied.
  Non-classical types require `calc_se = TRUE`; cluster types are
  unavailable for dynamic (4D) `W`.

- ...:

  Additional arguments (unused).

## Value

A matrix with one row per parameter and columns for the lower and upper
bounds, labeled by percentage (e.g., `"2.5 %"` and `"97.5 %"`).

## See also

[`boot_sir`](https://netify-dev.github.io/sir/reference/boot_sir.md) for
bootstrap inference,
[`vcov.sir`](https://netify-dev.github.io/sir/reference/vcov.sir.md) for
the variance-covariance matrix.
