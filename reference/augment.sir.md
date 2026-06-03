# Augment Data With SIR Model Fitted Values and Residuals

Returns the long-format dyad-time table of observed outcomes, fitted
values, and residuals from a fitted SIR model. Because SIR data are
arrays rather than a single data frame, the returned table is built from
the model's own `Y`/`fitted.values`/`residuals` with sender, receiver,
and time index columns.

## Usage

``` r
# S3 method for class 'sir'
augment(x, ...)
```

## Arguments

- x:

  A fitted `sir` object from
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md).

- ...:

  Unused, for generic compatibility.

## Value

A data frame with columns `sender`, `receiver`, `time`, `.observed`,
`.fitted`, `.resid` (response residual), and `.resid_pearson`. Diagonal
(self-tie) cells and any cells excluded from the likelihood are dropped.

## Examples

``` r
dat <- sim_sir(m = 10, T_len = 20, p = 2, q = 1, family = "poisson", seed = 1)
fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 1)
head(augment(fit))
#>   sender receiver time .observed   .fitted    .resid .resid_pearson
#> 2      2        1    1         1 2.0086935 -1.008694     -0.7117089
#> 3      3        1    1         1 0.8490400  0.150960      0.1638317
#> 4      4        1    1         4 1.4415275  2.558473      2.1309305
#> 5      5        1    1         4 0.9732095  3.026791      3.0681684
#> 6      6        1    1         2 0.7988776  1.201122      1.3438386
#> 7      7        1    1         3 1.5139874  1.486013      1.2077064
```
