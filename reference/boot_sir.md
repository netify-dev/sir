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
  type = c("block", "parametric", "dyad"),
  seed = NULL,
  trace = FALSE,
  cores = 1L
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

  Character. Inference type: `"block"` (default) resamples whole time
  periods as independent blocks, preserving within-period network
  dependence but not serial dependence; `"dyad"` is a delete-one-actor
  jackknife on the induced sub-network (captures dyadic dependence, and
  the only inference path for full-bilinear bipartite fits);
  `"parametric"` simulates new outcomes from the fitted model. `R` is
  ignored for `"dyad"` (it enumerates all actor deletions).

- seed:

  Optional integer for reproducibility. Sets the random seed before
  resampling. For exact reproducibility when `cores > 1`, the
  "L'Ecuyer-CMRG" RNG is used so parallel workers draw independent
  streams; results are then reproducible given the same `seed` and
  `cores`. Serial runs (`cores = 1`) with the same `seed` are always
  reproducible.

- trace:

  Logical. If TRUE, shows progress and prints a verbose message every 10
  serial replicates. Ignored when `cores > 1`.

- cores:

  Integer. Number of CPU cores to use. Default is 1 (serial). When
  greater than 1, replicates are run in parallel with
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html)
  (forking). Forking is unavailable on Windows, so there `cores > 1`
  falls back to serial with a one-time warning.

## Value

An object of class `"boot_sir"` with components:

- coefs:

  R x n_params matrix of bootstrap coefficient estimates. Rows for
  failed replicates contain NA.

- se:

  Named numeric vector of bootstrap standard errors (one per parameter).

- cov:

  For `type = "dyad"`, the jackknife variance-covariance matrix (NULL
  for the block/parametric bootstraps).

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

Three bootstrap strategies are available:

- block:

  Resamples whole time periods with replacement as independent period
  blocks. This preserves the within-period network dependence structure
  but not serial dependence across adjacent periods. Best when T is
  moderately large (T \>= 10).

- dyad:

  A delete-one-actor jackknife: each actor is dropped in turn and the
  model is refit on the induced sub-network, capturing the
  cross-sectional dyadic dependence (shared-sender / shared-receiver
  correlation) the block bootstrap ignores. Deletion (rather than
  resampling with replacement) avoids duplicating actors into the
  bilinear \\A X B'\\ sums, which would attenuate the influence
  coefficients toward zero. For a square directed network the same actor
  is dropped from both axes together; for a bipartite network senders
  and receivers are dropped separately and the two one-way jackknife
  covariances are summed. Standard errors come from the jackknife
  covariance and intervals are normal (`estimate +/- z * se`). This is
  also the estimator
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md) reuses
  automatically when `calc_se = TRUE` cannot form analytic SEs (a
  singular/ill-conditioned Hessian or a full-bilinear bipartite fit): it
  attaches the same `$cov` so `vcov`/`confint`/`tidy` return jackknife
  inference without an explicit `boot_sir` call.

- parametric:

  Simulates new outcome arrays from the fitted model using the estimated
  parameters and the specified family distribution. Better when T is
  small but the model is well-specified.

For `block` and `parametric`, each replicate refits the full SIR model;
replicates that fail to converge are dropped and reported, standard
errors are the column standard deviations of the successful replicates,
and confidence intervals use the percentile method. The `dyad` jackknife
instead enumerates all delete-one-actor refits and reports the jackknife
standard error and a normal interval.

## See also

[`confint.sir`](https://netify-dev.github.io/sir/reference/confint.sir.md)
to use bootstrap intervals,
[`confint.boot_sir`](https://netify-dev.github.io/sir/reference/confint.boot_sir.md)
for direct interval extraction.

## Examples

``` r
# \donttest{
dat <- sim_sir(m = 8, T_len = 12, p = 2, q = 1, family = "poisson", seed = 1)
model <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z,
             family = "poisson", seed = 1)

# use a larger R (for example 500+) for publication-quality intervals
boot_result <- boot_sir(model, R = 10, seed = 42)
print(boot_result)
#> 
#> Bootstrap SIR Results
#> Type: block | Replicates: 10/10 valid
#> 
#>             Estimate Boot SE   2.5 %  97.5 %
#> (Z) Z1        0.3885  0.0699  0.2591  0.4728
#> (alphaW) W2  -0.0452  0.1685 -0.3325  0.1841
#> (betaW) W1    0.0384  0.0863 -0.0797  0.1727
#> (betaW) W2   -0.3911  0.0471 -0.4335 -0.2918
#> 

# use bootstrap CIs with confint
confint(model, boot = boot_result)
#>                   2.5 %     97.5 %
#> (Z) Z1       0.25912539  0.4727739
#> (alphaW) W2 -0.33252941  0.1840802
#> (betaW) W1  -0.07972307  0.1727159
#> (betaW) W2  -0.43347988 -0.2917629
# }
```
