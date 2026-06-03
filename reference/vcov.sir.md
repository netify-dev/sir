# Variance-Covariance Matrix for SIR Model Parameters

Returns the variance-covariance matrix of the estimated parameters.
Several types are available; for relational (dyadic) data the
cluster-robust types are usually more cautious because the classical and
HC0-robust covariances assume independent dyad-periods and often
undercover.

## Usage

``` r
# S3 method for class 'sir'
vcov(object, type = c("cluster", "classical", "robust", "twoway", "dyad"), ...)
```

## Arguments

- object:

  A fitted `sir` object from
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md).

- type:

  Character string:

  - `"cluster"` (default) — multiway cluster-robust covariance on
    sender, receiver, and time margins for directed-network data when
    the Hessian bread is stable.

  - `"classical"` — inverse-Hessian covariance.

  - `"robust"` — HC0 sandwich; corrects heteroskedasticity /
    overdispersion only, *not* dyadic dependence.

  - `"twoway"` — alias for `"cluster"`.

  - `"dyad"` — clusters the directed dyad across time.

  The cluster types require the classical covariance as their bread
  (`calc_se = TRUE`, the default) and are unavailable for dynamic (4D)
  `W` and for full-bilinear bipartite fits (use
  `boot_sir(type = "dyad")` there).

- ...:

  Additional arguments (unused).

## Value

A square matrix with rows and columns named by parameter. Returns NULL
if standard errors were not computed (`calc_se = FALSE`).

## See also

[`confint.sir`](https://netify-dev.github.io/sir/reference/confint.sir.md)
for confidence intervals,
[`boot_sir`](https://netify-dev.github.io/sir/reference/boot_sir.md) for
resampling-based inference.
