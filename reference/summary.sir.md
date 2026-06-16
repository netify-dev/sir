# Summary of a Fitted SIR Model

Produces a detailed summary of the fitted model including coefficient
estimates with standard errors and p-values, model fit statistics
(log-likelihood, AIC, BIC), convergence status, and summaries of the
estimated influence matrices A and B.

## Usage

``` r
# S3 method for class 'sir'
summary(object, ...)
```

## Arguments

- object:

  A fitted `sir` object from
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md).

- ...:

  Additional arguments (unused).

## Value

An object of class `"summary.sir"` containing:

- coefficients:

  Data frame with columns `coef`, `se` (classical SE), `rse`
  (cluster-robust SE), `t_se` (classical z, the printed "z value"),
  `t_rse` (robust z), `p.value`, and `sig` (significance codes). The
  printed z/p/stars are classical.

- loglik:

  Log-likelihood at convergence.

- aic:

  AIC value.

- bic:

  BIC value.

- converged:

  Logical convergence indicator.

- iterations:

  Iteration count.

- A.summary:

  List with mean, sd, and range of off-diagonal entries in the sender
  effects matrix.

- B.summary:

  Same for the receiver effects matrix.

## See also

[`print.summary.sir`](https://netify-dev.github.io/sir/reference/print.summary.sir.md)
for the printed output.
