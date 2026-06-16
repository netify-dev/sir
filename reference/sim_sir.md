# Simulate Data from a Social Influence Regression Model

Generates synthetic network data from a known SIR data-generating
process. Useful for testing, benchmarking, and pedagogical
demonstrations. The function simulates Y\[i,j,t\] from the specified
family using influence covariates W, lagged network state X, and
optional exogenous covariates Z.

## Usage

``` r
sim_sir(
  m,
  T_len,
  p = 2,
  q = 1,
  family = "poisson",
  alpha = NULL,
  beta = NULL,
  theta = NULL,
  W = NULL,
  sigma = 1,
  seed = NULL,
  symmetric = FALSE,
  gain = NULL,
  ...
)
```

## Arguments

- m:

  Integer. Number of nodes in the network.

- T_len:

  Integer. Number of time periods.

- p:

  Integer. Number of influence covariates in W. Default is 2.

- q:

  Integer. Number of exogenous covariates in Z. Default is 1. Set to 0
  for no exogenous covariates.

- family:

  Character string: `"poisson"` (default), `"normal"`, or `"binomial"`.

- alpha:

  Numeric vector of length p for sender influence weights. For a
  directed process the first element (alpha_1) is fixed at 1 for
  identifiability; for `symmetric = TRUE` all elements are kept as given
  (the anchor need not be 1). If NULL (default), drawn from N(0, 0.3)
  with alpha_1 = 1. Use `seed` for reproducibility.

- beta:

  Numeric vector of length p for receiver influence weights. If NULL
  (default), drawn from N(0, 0.3).

- theta:

  Numeric vector of length q for exogenous covariate effects. If NULL
  (default), drawn from U(-0.5, 0.5).

- W:

  Optional 3D array (m x m x p) of influence covariates. If NULL
  (default), generated with standard normal entries (a dense,
  well-conditioned influence design). The diagonal is set to zero
  because self-ties are not part of the one-mode SIR likelihood.

- sigma:

  Numeric. Standard deviation for the normal family. Default 1.

- seed:

  Optional integer for reproducibility. When supplied, the seed is set
  locally and the caller's global RNG state is restored on exit, so a
  subsequent draw (e.g. `runif`) in the caller is left unperturbed.

- symmetric:

  Logical. If TRUE, simulate a genuine **undirected** network: `W` is
  symmetrized, `beta` is tied to `alpha` so `B = A`, and each `Y`/`Z`
  slice is symmetric. For the count/continuous recursion the influence
  covariates are rescaled (when auto-generated) so the symmetric gain
  \\\rho(A)^2/(m-1)\\ stays below 1. Pairs with
  `sir(..., symmetric = TRUE)`, which fits and reports the shared
  operator as `gamma`; `sim_sir` stores that same vector in `$alpha`
  (which equals `$beta` here). Default FALSE.

- gain:

  Optional numeric in (0, 1). Target spectral gain
  \\\rho(A)\rho(B)/(m-1)\\ for the Poisson/Normal lagged recursion. When
  supplied, the influence operators are rescaled to hit it exactly
  (scaling `beta`/`B` for directed fits, the shared operator for
  symmetric, leaving `alpha_1 = 1` intact), so the influence term
  carries a chosen, strong-but-stationary share of the dynamics. Use a
  larger value (e.g. 0.9) for a strong, recoverable influence signal;
  values near 1 approach non-stationarity. NULL (default) keeps the
  conservative auto-rescaling. Ignored for `family = "binomial"`.

- ...:

  Unused; catches mistyped arguments and reports a clear error.

## Value

A list with components:

- Y:

  3D array (m x m x T_len) of simulated outcomes.

- W:

  3D array (m x m x p) of influence covariates.

- X:

  3D array (m x m x T_len) of the (scaled) lagged network state used in
  the bilinear mean: `X[,,t]` is `log(Y[,,t-1] + 1)` (Poisson) or
  `Y[,,t-1]` (otherwise), divided by `(m - 1)`.

- Z:

  4D array (m x m x q x T_len) of exogenous covariates, or NULL if q =
  0.

- alpha:

  True alpha vector (length p; alpha_1 = 1 for directed, kept as
  supplied for symmetric). For symmetric fits this equals `beta` and is
  the shared operator reported by the fit as gamma.

- beta:

  True beta vector (length p).

- theta:

  True theta vector (length q).

- A:

  True sender influence matrix (m x m).

- B:

  True receiver influence matrix (m x m).

- family:

  The distribution family used.

## Details

The influence-carrying state `X` is scaled by `1 / (m - 1)` before it
enters the bilinear mean. The bilinear term \\A X_t B^\top\\ sums over
all \\(m-1)\\ off-diagonal partners on each side, so without this
scaling the linear predictor would grow with network size and (for the
Poisson/log link) the conditional mean would saturate, making the
influence parameters unrecoverable. The returned `X` already includes
this scaling, so a plain
[`sir`](https://netify-dev.github.io/sir/reference/sir.md) fit on
`(Y, W, X, Z)` recovers the same `A`, `B` used to generate the data.

For the Poisson and Normal families the lagged recursion is stationary
only when the spectral gain \\\rho(A)\rho(B)/(m-1)\\ is below 1. By
default, auto-generated coefficients are rescaled to a conservative gain
(at most 0.8) so simulations are reliably stable, and user-supplied
coefficients are used exactly as given (a gain \\\ge 1\\ triggers an
explosive-series warning). Set `gain` to target a specific value: the
influence operators are then rescaled so \\\rho(A)\rho(B)/(m-1)\\ equals
`gain` exactly, letting the bilinear term carry a chosen,
strong-but-stationary share of the dynamics. A larger `gain` (say 0.9)
makes the influence mechanism dominate the lagged dynamics while still
recovering cleanly; values near 1 approach non-stationarity. Binomial
outcomes are bounded, so `gain` does not apply.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate Poisson network and recover the parameters
dat <- sim_sir(m = 15, T_len = 30, p = 2, q = 1, family = "poisson", seed = 42)
fit <- sir(dat$Y, dat$W, dat$X, dat$Z, family = "poisson")
cbind(true = c(dat$theta, dat$alpha[-1], dat$beta), estimated = coef(fit))
} # }
```
