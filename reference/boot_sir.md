# Bootstrap Inference for SIR Model Parameters

Computes bootstrap standard errors and confidence intervals for SIR
model parameters. This is the recommended approach for inference when
the Hessian is singular or ill-conditioned, which is common in models
with bilinear influence terms (i.e., when `fix_receiver = FALSE`).

## Usage

``` r
boot_sir(
  sir_fit,
  R = 200,
  type = c("block", "parametric"),
  seed = NULL,
  trace = FALSE
)
```

## Arguments

- sir_fit:

  A fitted `sir` object from
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md).

- R:

  Integer. Number of bootstrap replicates. Default is 200. Increase to
  500-1000 for publication-quality intervals.

- type:

  Character. Bootstrap type: `"block"` (default) resamples time periods
  with replacement; `"parametric"` simulates new outcomes from the
  fitted model.

- seed:

  Optional integer for reproducibility. Sets the random seed before
  resampling.

- trace:

  Logical. If TRUE, prints progress every 10 replicates.

## Value

An object of class `"boot_sir"` with components:

- coefs:

  R x n_params matrix of bootstrap coefficient estimates. Rows for
  failed replicates contain NA.

- se:

  Named numeric vector of bootstrap standard errors (one per parameter).

- ci_lo:

  Lower 2.5% percentile bounds.

- ci_hi:

  Upper 97.5% percentile bounds.

- point_est:

  Point estimates from the original fit.

- param_names:

  Character vector of parameter names.

- n_valid:

  Number of successful bootstrap replicates.

- n_total:

  Total number of replicates attempted.

- type:

  The bootstrap type used.

- family:

  The distribution family.

## Details

Two bootstrap strategies are available:

- block:

  Resamples time periods with replacement. Preserves the within-period
  dependence structure. Best when T is moderately large (T \>= 10).

- parametric:

  Simulates new outcome arrays from the fitted model using the estimated
  parameters and the specified family distribution. Better when T is
  small but the model is well-specified.

Each replicate refits the full SIR model. Replicates that fail to
converge are dropped and reported. Standard errors are column standard
deviations of the successful replicates. Confidence intervals use the
percentile method.

## See also

[`confint.sir`](https://netify-dev.github.io/sir/reference/confint.sir.md)
to use bootstrap intervals,
[`confint.boot_sir`](https://netify-dev.github.io/sir/reference/confint.boot_sir.md)
for direct interval extraction.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- sir(Y, W, X, family = "poisson")

# block bootstrap with 200 replicates
boot_result <- boot_sir(model, R = 200, seed = 42)
print(boot_result)

# use bootstrap CIs with confint
confint(model, boot = boot_result)

# parametric bootstrap
boot_par <- boot_sir(model, R = 200, type = "parametric")
} # }
```
