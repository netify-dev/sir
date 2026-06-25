# Variance-Covariance Matrix for SIR Model Parameters

Returns the variance-covariance matrix of the estimated parameters.
Several types are available; for relational (dyadic) data the
cluster-robust types are usually more cautious because the classical and
HC0-robust covariances assume independent dyad-periods and often
undercover.

## Usage

``` r
# S3 method for class 'sir'
vcov(object, type = c("cluster", "classical", "robust"), ...)
```

## Arguments

- object:

  A fitted `sir` object from
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md).

- type:

  Character string:

  - `"cluster"` (default) — actor-clustered cluster-robust covariance:
    each cell's score is stacked onto both endpoint actors and summed
    within actor, with an HC1 small-sample factor. This is a
    conservative actor-margin sandwich (not textbook two-way CGM), used
    identically for directed and symmetric fits. The returned matrix
    carries a `"cluster_df"` attribute (number of actor clusters minus
    one) that `confint`/`tidy` use for the `t(G - 1)` reference.

  - `"classical"` — inverse-Hessian covariance; valid when dyads are
    independent, and tighter, but it undercovers under dyadic
    dependence.

  - `"robust"` — the HC0 sandwich for directed fits; for symmetric fits
    (which have no separate HC0 path) it aliases to `"cluster"`.

  The cluster type requires the classical covariance as its bread
  (`calc_se = TRUE`, the default) and is supported for directed,
  symmetric, and dynamic (4D) `W` fits; it is unavailable only for
  full-bilinear bipartite fits. For that case (and any fit whose Hessian
  was ill-conditioned), `sir` attaches a delete-one-actor jackknife
  covariance automatically: when `object$se_source == "jackknife"`,
  `vcov` returns that jackknife covariance for *every* `type` (the
  `type` argument is ignored, since no analytic classical/HC0/cluster
  covariance exists). `confint` pairs the cluster type with a `t(G - 1)`
  reference (G = number of actors).

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
