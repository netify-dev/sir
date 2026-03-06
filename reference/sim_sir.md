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
  seed = NULL
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
  are free. If NULL (default), drawn from N(0, 0.3).

- beta:

  Numeric vector of length p for receiver influence weights. If NULL
  (default), drawn from N(0, 0.3).

- theta:

  Numeric vector of length q for exogenous covariate effects. If NULL
  (default), drawn from U(-0.5, 0.5).

- W:

  Optional 3D array (m x m x p) of influence covariates. If NULL
  (default), generated with standard normal entries.

- sigma:

  Numeric. Standard deviation for the normal family. Default 1.

- seed:

  Optional integer for reproducibility.

## Value

A list with components:

- Y:

  3D array (m x m x T_len) of simulated outcomes.

- W:

  3D array (m x m x p) of influence covariates.

- X:

  3D array (m x m x T_len) of lagged network state.

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

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate Poisson network
dat <- sim_sir(m = 15, T_len = 10, p = 2, q = 1, family = "poisson", seed = 42)

# Fit model to recover parameters
fit <- sir(dat$Y, dat$W, dat$X, dat$Z, family = "poisson")
cbind(true = c(dat$theta, dat$alpha[-1], dat$beta), estimated = coef(fit))
} # }
```
