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

  Numeric vector of length p for sender influence weights. The first
  element (alpha_1) is fixed at 1 for identifiability; only alpha_2:p
  are free. If NULL (default), drawn from N(0, 0.3). Use `seed` for
  reproducibility.

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

  True alpha vector (length p, with alpha_1 = 1).

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

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate Poisson network and recover the parameters
dat <- sim_sir(m = 15, T_len = 30, p = 2, q = 1, family = "poisson", seed = 42)
fit <- sir(dat$Y, dat$W, dat$X, dat$Z, family = "poisson")
cbind(true = c(dat$theta, dat$alpha[-1], dat$beta), estimated = coef(fit))
} # }
```
