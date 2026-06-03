# One-Row Model Summary for a SIR Model

Returns a single-row data frame of model-level statistics, in the layout
broom consumers expect for model comparison tables.

## Usage

``` r
# S3 method for class 'sir'
glance(x, ...)
```

## Arguments

- x:

  A fitted `sir` object from
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md).

- ...:

  Unused, for generic compatibility.

## Value

A one-row data frame with columns `nobs`, `df` (number of estimated
parameters), `logLik`, `AIC`, `BIC`, `n_nodes`, `n_periods`,
`n_influence_covar` (p), `family`, `method`, and `converged`.

## Examples

``` r
dat <- sim_sir(m = 10, T_len = 20, p = 2, q = 1, family = "poisson", seed = 1)
fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 1)
glance(fit)
#>   nobs df    logLik      AIC      BIC n_nodes n_periods n_influence_covar
#> 1 1800  4 -2412.457 4832.913 4854.895      10        20                 2
#>    family method converged
#> 1 poisson    ALS      TRUE
```
