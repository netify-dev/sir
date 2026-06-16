# Score predicted networks against observed outcomes

Computes out-of-sample scores comparing predicted to actual networks. By
default, square arrays are scored on the off-diagonal non-missing cells
because one-mode self-ties are not modeled; set `drop_diagonal = FALSE`
for square bipartite arrays whose diagonal cells are real
sender-receiver observations.

## Usage

``` r
score_sir(
  actual,
  predicted,
  family = c("poisson", "normal", "binomial"),
  drop_diagonal = TRUE
)
```

## Arguments

- actual, predicted:

  numeric arrays or matrices of identical shape, on the response scale
  (counts, probabilities, or means).

- family:

  character; one of `"poisson"`, `"normal"`, `"binomial"`.

- drop_diagonal:

  logical; if TRUE (default), exclude diagonal cells for square
  matrices/arrays before scoring.

## Value

a named numeric vector of scores. `normal`: rmse, mae. `poisson`: rmse,
mae, deviance. `binomial`: rmse, mae, logloss, brier, auc. Lower is
better for every score except `auc` (higher is better). Notes:
`deviance` is the per-cell *mean* unit deviance (not the total, so it is
comparable across origins of different size); non-finite `actual` cells
are excluded, while non-finite `predicted` values on scored cells raise
an error; `logloss` clamps predicted probabilities away from 0/1, so
against a hard-label (0/1) baseline it is dominated by that clamp –
prefer `brier`/`auc` there; `auc` is `NA` when the held-out cells
contain a single class. Binary `actual` must be coded 0/1 and poisson
`actual` must be non-negative integer counts, else an error is raised.

## See also

[`cv_sir`](https://netify-dev.github.io/sir/reference/cv_sir.md) for
cross-validated scoring,
[`forecast.sir_fit`](https://netify-dev.github.io/sir/reference/forecast.sir_fit.md)
for producing the predictions to score.

## Examples

``` r
a <- matrix(rpois(100, 2), 10, 10); diag(a) <- NA
p <- a + 0.5
score_sir(a, p, "poisson")
#>      rmse       mae  deviance 
#> 0.5000000 0.5000000 0.1972497 
```
