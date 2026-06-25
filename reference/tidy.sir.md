# Tidy a SIR Model into a Data Frame of Coefficients

Returns a data frame with one row per estimated parameter, in the layout
expected by broom consumers such as modelsummary and gtsummary. Each
parameter is tagged by its role (`component`: exogenous `theta`, sender
`alpha`, receiver `beta`, or shared `gamma` for a symmetric fit).

## Usage

``` r
# S3 method for class 'sir'
tidy(
  x,
  conf.int = FALSE,
  conf.level = 0.95,
  se.type = c("cluster", "classical", "robust"),
  ...
)
```

## Arguments

- x:

  A fitted `sir` object from
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md).

- conf.int:

  Logical; if `TRUE`, add `conf.low`/`conf.high` Wald interval columns.
  Default `FALSE`.

- conf.level:

  Confidence level for the interval. Default 0.95.

- se.type:

  Which standard errors to report: `"cluster"` (default; actor-clustered
  sandwich, each cell scored onto both endpoint actors, for directed,
  symmetric, and dynamic (4D) `W` fits), `"classical"`
  (inverse-Hessian), or `"robust"` (HC0 sandwich). Ignored when the fit
  carries `se_source == "jackknife"` (analytic SEs were unavailable, so
  the delete-one-actor jackknife standard errors are reported for every
  type).

- ...:

  Unused, for generic compatibility.

## Value

A data frame with columns `term`, `component`, `estimate`, `std.error`,
`statistic`, `p.value`, and (optionally) `conf.low`, `conf.high`.

## Examples

``` r
dat <- sim_sir(m = 10, T_len = 20, p = 2, q = 1, family = "poisson", seed = 1)
fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 1)
tidy(fit)
#>          term component   estimate  std.error  statistic      p.value
#> 1      (Z) Z1     theta  0.4478604 0.02438316 18.3676096 1.921382e-08
#> 2 (alphaW) W2     alpha -0.1117818 0.28001836 -0.3991945 6.990566e-01
#> 3  (betaW) W1      beta  0.1044978 0.07596178  1.3756629 2.021912e-01
#> 4  (betaW) W2      beta -0.1958658 0.04609193 -4.2494595 2.144178e-03
tidy(fit, conf.int = TRUE)
#>          term component   estimate  std.error  statistic      p.value
#> 1      (Z) Z1     theta  0.4478604 0.02438316 18.3676096 1.921382e-08
#> 2 (alphaW) W2     alpha -0.1117818 0.28001836 -0.3991945 6.990566e-01
#> 3  (betaW) W1      beta  0.1044978 0.07596178  1.3756629 2.021912e-01
#> 4  (betaW) W2      beta -0.1958658 0.04609193 -4.2494595 2.144178e-03
#>      conf.low  conf.high
#> 1  0.39270186  0.5030190
#> 2 -0.74522733  0.5216637
#> 3 -0.06733968  0.2763353
#> 4 -0.30013298 -0.0915986
```
