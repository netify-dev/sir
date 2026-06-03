# Tidy a SIR Model into a Data Frame of Coefficients

Returns a data frame with one row per estimated parameter, in the layout
expected by broom consumers such as modelsummary and gtsummary. Each
parameter is tagged by its role (`component`: exogenous `theta`, sender
`alpha`, or receiver `beta`).

## Usage

``` r
# S3 method for class 'sir'
tidy(
  x,
  conf.int = FALSE,
  conf.level = 0.95,
  se.type = c("cluster", "classical", "robust", "twoway", "dyad"),
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

  Which standard errors to report: `"cluster"` (default; multiway
  sender, receiver, and time clustering for supported static fits),
  `"classical"` (inverse-Hessian), `"robust"` (HC0 sandwich), `"twoway"`
  (alias for `"cluster"`), or `"dyad"` (directed-dyad clustering).

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
#>          term component   estimate   std.error  statistic      p.value
#> 1      (Z) Z1     theta  0.4478604 0.009540935 46.9409338 0.000000e+00
#> 2 (alphaW) W2     alpha -0.1117818 0.220142366 -0.5077705 6.116143e-01
#> 3  (betaW) W1      beta  0.1044978 0.052961170  1.9731022 4.848392e-02
#> 4  (betaW) W2      beta -0.1958658 0.021712466 -9.0208913 1.865651e-19
tidy(fit, conf.int = TRUE)
#>          term component   estimate   std.error  statistic      p.value
#> 1      (Z) Z1     theta  0.4478604 0.009540935 46.9409338 0.000000e+00
#> 2 (alphaW) W2     alpha -0.1117818 0.220142366 -0.5077705 6.116143e-01
#> 3  (betaW) W1      beta  0.1044978 0.052961170  1.9731022 4.848392e-02
#> 4  (betaW) W2      beta -0.1958658 0.021712466 -9.0208913 1.865651e-19
#>        conf.low  conf.high
#> 1  0.4291605186  0.4665603
#> 2 -0.5432529100  0.3196893
#> 3  0.0006958146  0.2082998
#> 4 -0.2384214428 -0.1533101
```
