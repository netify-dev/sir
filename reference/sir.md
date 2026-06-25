# Social Influence Regression (SIR) Model

Fits a Social Influence Regression model for network data with lagged,
model-implied association channels. The SIR model captures how prior
network states predict later outcomes through bilinear interaction
terms, allowing for both sender and receiver channels in directed
networks. Causal interpretation requires additional research-design
assumptions.

The model decomposes network influence into two components:

- **Sender influence (A matrix)**: A\[i,k\] measures how much node k's
  behavior (via X) shapes node i's outgoing ties.

- **Receiver influence (B matrix)**: B\[j,l\] measures how node l's
  position shapes node j's incoming ties.

## Usage

``` r
sir(
  Y,
  W = NULL,
  X = NULL,
  Z = NULL,
  family,
  method = "ALS",
  calc_se = TRUE,
  fix_receiver = FALSE,
  symmetric = FALSE,
  bipartite = NULL,
  W_recv = NULL,
  kron_mode = FALSE,
  seed = NULL,
  ...
)
```

## Arguments

- Y:

  A three-dimensional array of dimensions (m x m x T) containing the
  network outcomes. Y\[i,j,t\] represents the directed outcome from node
  i to node j at time t. Can contain NA values for missing observations.
  For one-mode square networks, self-tie diagonals are excluded from
  fitting, likelihood evaluation, and reported observation counts.

- W:

  Optional influence covariate array, either:

  - **3D array** (m x m x p): Static influence covariates. W\[i,k,r\]
    represents the r-th covariate for the channel by which source node k
    can shape target node i; the same orientation is used for the
    receiver-side matrix B\[j,l\]. Directed W slices therefore need the
    same target-by-source orientation as A and B, not necessarily the
    same sender-to-receiver orientation as Y\[i,j,t\]. The same W is
    used for all time periods.

  - **4D array** (m x m x p x T): Dynamic (time-varying) influence
    covariates. W\[i,k,r,t\] allows the influence structure to change
    over time. Parameters (alpha, beta) are still estimated jointly
    across all periods, but the influence matrices A_t and B_t vary
    with t. Only ALS method is supported for 4D W.

  Common choices include graph Laplacians, geographic distance matrices,
  or node-level covariates expanded to edge-level. For one-mode square
  networks, W diagonals are set to zero so self-influence channels are
  not estimated. If NULL or p=0, no network influence structure is
  included.

- X:

  Optional three-dimensional array of dimensions (m x m x T)
  representing the network state that carries influence. Typically this
  is a lagged version of Y (e.g., X\[,,t\] = Y\[,,t-1\]). If NULL and W
  is provided, an error is thrown. X determines which network patterns
  influence future outcomes.

- Z:

  Optional array of exogenous covariates. Can be either:

  - 3D array (m x m x T): Single covariate varying across edges and time

  - 4D array (m x m x q x T): Multiple (q) covariates

  Examples include dyadic covariates (trade agreements, geographic
  distance) or node-level attributes (GDP, population) expanded to
  edge-level.

- family:

  Character string specifying the distribution family and link function.
  Must be one of "poisson", "normal", or "binomial". The choice depends
  on the nature of your outcome variable.

- method:

  Character string specifying the estimation method. Either `"ALS"` (the
  default alternating GLM/IRLS engine; the name is retained for API
  compatibility) or `"optim"` (direct optimization via BFGS).

- calc_se:

  Logical indicating whether to calculate standard errors for the
  parameters (default TRUE). Standard errors are computed from the
  observed information matrix. If those analytic standard errors cannot
  be formed — a singular or ill-conditioned Hessian, or a path with no
  closed-form covariance such as full-bilinear bipartite fits — and the
  model converged, `sir` automatically falls back to the
  delete-one-actor jackknife covariance (reproducible for a given
  `seed`), so [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`summary()`](https://rdrr.io/r/base/summary.html), and
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) still
  return standard errors. The fit then carries `se_source = "jackknife"`
  and prints a one-line note. Set FALSE to skip standard errors entirely
  (and the fallback) when they are not needed.

- fix_receiver:

  Logical. If TRUE, fixes B = I (identity matrix) and estimates only
  (theta, alpha). This eliminates the bilinear identification problem
  (scaling ambiguity between A and B) by removing the receiver influence
  channel. The model becomes a standard GLM, yielding model-based
  standard errors under the usual GLM assumptions. Appropriate when
  receiver effects are negligible. Default is FALSE.

- symmetric:

  Logical. If TRUE, fits the genuine **undirected** model with a single
  shared influence operator \\A = B = \sum_k \gamma_k W_k\\, so the
  bilinear term is the quadratic form \\A X A'\\, symmetric in \\(i,j)\\
  by construction. All \\\gamma_k\\ are estimated (the quadratic form's
  scale is identified by the data); \\A\\ is identified only up to its
  overall sign (\\A X A' = (-A) X (-A)'\\), fixed so the
  largest-magnitude \\\gamma_k\\ is positive. Requires a square network
  and symmetric, zero-diagonal influence covariates `W` (non-zero W
  diagonals are zeroed, since the quadratic form must be
  self-feedback-free); static (3D) or dynamic (4D, time-varying) `W` are
  both supported, the dynamic case giving a per-period operator \\A_t
  X_t A_t'\\. Continuous asymmetric `Y` is averaged across triangles;
  Poisson/Binomial `Y` must already be symmetric. Estimation is BFGS
  with an analytic gradient over the upper-triangle off-diagonal cells;
  [`confint()`](https://rdrr.io/r/stats/confint.html)/[`vcov()`](https://rdrr.io/r/stats/vcov.html)/
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) default to
  the actor-clustered cluster-robust SE (the same estimator used for
  directed fits, with a \\t(G-1)\\ reference), while
  [`summary()`](https://rdrr.io/r/base/summary.html) prints the
  classical SE. Request classical intervals with `se.type = "classical"`
  (see
  [`confint.sir`](https://netify-dev.github.io/sir/reference/confint.sir.md)).
  Cannot be combined with `fix_receiver`. Default is FALSE.

- bipartite:

  Logical or NULL. Indicates whether the network is bipartite (senders
  and receivers are distinct node sets). If NULL (the default),
  bipartite status is inferred from Y: non-square arrays (n1 != n2) are
  treated as bipartite. Set to TRUE explicitly for square arrays where
  senders and receivers are nonetheless distinct populations. Setting
  FALSE on a non-square Y raises an error. Bipartite networks require
  `fix_receiver = TRUE` unless a separate `W_recv` is supplied (see
  below).

- W_recv:

  Optional receiver-side influence covariate array (n2 x n2 x p2) for a
  **full-bilinear bipartite** fit. When supplied, the sender influence
  `A` is built from `W` (n1 x n1 x p) and the receiver influence `B`
  from `W_recv`, fitting the complete \\A X B'\\ model for two-mode data
  via alternating GLM (rather than collapsing to `B = I`). `alpha_1 = 1`
  pins the scale; all `beta` are free. Works for both rectangular (n1 !=
  n2) and square two-mode networks. This path uses a pure-R
  alternating-GLM estimator (it does not call the C++ likelihood
  kernel), so it is slower than the square one-mode path for large
  networks; it draws random restarts (see `n_restarts` in `...`,
  default 5) and respects `seed`. Analytic standard errors are not
  available; use
  [`boot_sir`](https://netify-dev.github.io/sir/reference/boot_sir.md)
  with `type = "dyad"`. Default NULL.

- kron_mode:

  Logical. **Not yet implemented** (reserved for a future release):
  would estimate an unconstrained p x p coefficient matrix C instead of
  the rank-at-most-one alpha beta' factorization. Setting
  `kron_mode = TRUE` currently raises an error. Default is FALSE.

- seed:

  Optional integer. If supplied, sets the random seed before fitting so
  that runs are reproducible (the estimators use random starting
  values). The global RNG state is restored on exit. Default NULL.

- ...:

  Additional arguments passed to the fitting functions:

  - `trace`: Logical or integer controlling output verbosity.

  - `tol`: Convergence tolerance for ALS (default 1e-8).

  - `max_iter`: Maximum ALS iterations (default 100).

  - `n_restarts`: Random restarts for the full-bilinear bipartite
    (`W_recv`) estimator (default 5); the lowest-deviance fit is kept.

## Value

An object of class `"sir"` with the following components:

- summ:

  Data frame of parameter estimates with columns `coef`, `se` (classical
  SE), `rse` (robust/sandwich SE), `t_se` (z-statistic using classical
  SE), `t_rse` (z-statistic using robust SE). Row names identify each
  parameter. The `rse`/`t_rse` columns are `NA` for fits with no
  separate HC0 path (notably symmetric/undirected fits; all SE columns
  are `NA` for full-bilinear bipartite fits). The default cluster-robust
  SEs/intervals come from `sqrt(diag(vcov(fit)))` and `confint(fit)`,
  not from `summ`.

- A:

  Sender influence matrix. For static W: n1 x n1 matrix. For dynamic
  (4D) W: n1 x n1 x T array. Off-diagonal entry A\[i,k\] measures how
  much node k's behavior (via X) shapes node i's outgoing ties. Diagonal
  is set to zero for one-mode square fits; for a full-bilinear bipartite
  fit (`W_recv`) A is the dense `sum_k alpha_k W[,,k]` with all
  sender-side entries retained.

- B:

  Receiver influence matrix. Identity when `fix_receiver = TRUE`;
  `n1 x n1` (same shape as A) for the square model; `n2 x n2` built from
  `W_recv` for a full-bilinear bipartite fit.

- tab:

  Numeric vector of all estimated parameters in order: \[theta_1, ...,
  theta_q, alpha_2, ..., alpha_p, beta_1, ..., beta_p\]. When
  `fix_receiver = TRUE`: \[theta_1, ..., theta_q, alpha_1, ...,
  alpha_p\]. For a symmetric fit: \[theta_1, ..., theta_q, gamma_1, ...,
  gamma_p\] (all gamma estimated; A is identified up to global sign,
  fixed so the largest-magnitude gamma is positive). For a full-bilinear
  bipartite fit: \[theta_1, ..., theta_q, alpha_2, ..., alpha_p, beta_1,
  ..., beta_p2\].

- theta:

  Coefficients for exogenous covariates Z (length q).

- alpha:

  Full alpha vector including the fixed alpha_1 = 1 (length p). When
  `fix_receiver = TRUE`, all alpha are free.

- beta:

  Coefficients for receiver influence covariates (length p, or length p2
  for a full-bilinear bipartite fit from `W_recv`). Empty when
  `fix_receiver = TRUE`.

- p2:

  Number of receiver-side influence covariates (only present for a
  full-bilinear bipartite fit).

- ll:

  Log-likelihood at convergence.

- family:

  The distribution family used (`"poisson"`, `"normal"`, or
  `"binomial"`).

- method:

  The estimation method used (`"ALS"` or `"optim"`).

- p:

  Number of influence covariates in W.

- q:

  Number of exogenous covariates in Z.

- m:

  Number of sender nodes (same as n1).

- n1:

  Number of sender (row) nodes.

- n2:

  Number of receiver (column) nodes.

- bipartite:

  Logical, TRUE if the network is treated as bipartite (non-square `Y`,
  `bipartite = TRUE`, or a `W_recv` fit).

- full_bilinear:

  Logical, TRUE for a full-bilinear bipartite fit (`W_recv` supplied);
  such fits have no analytic SEs.

- n_periods:

  Number of time periods.

- nobs:

  Number of non-missing observations used in estimation.

- fitted.values:

  Array (n1 x n2 x T) of fitted values on the response scale (counts for
  Poisson, probabilities for binomial, means for normal).

- residuals:

  List with three components: `response` (Y - fitted), `pearson`
  (standardized by variance function), and `deviance` (signed square
  root of deviance contributions).

- vcov:

  Variance-covariance matrix of parameters from the Hessian (classical
  SEs). NULL if `calc_se = FALSE`.

- vcov_robust:

  Sandwich (robust) variance-covariance matrix. NULL if
  `calc_se = FALSE` or computation failed.

- Y:

  The outcome array as used in fitting (with NAs from symmetric masking
  or Z missingness applied).

- W:

  The influence covariate array.

- X:

  The network state array (NAs replaced with 0).

- Z:

  The exogenous covariate array (converted to 4D if 3D input).

- fix_receiver:

  Logical, whether receiver effects were fixed.

- symmetric:

  Logical, whether the network was treated as undirected.

- kron_mode:

  Logical, whether Kronecker mode was used.

- iterations:

  Number of iterations until convergence.

- history:

  List with matrices ALPHA, BETA, THETA, DEV tracking parameter
  trajectories across iterations (useful for convergence diagnostics).

- convergence:

  Logical, TRUE if the algorithm converged.

- call:

  The matched function call.

- sigma2:

  Estimated error variance (only for `family = "normal"`).

- se_reliable:

  Logical, FALSE if the Hessian was ill-conditioned so the classical SEs
  should be treated with caution.

- se_source:

  Character flag for the reported standard errors: `"jackknife"` when
  analytic SEs could not be formed and `sir` fell back to the
  delete-one-actor jackknife (see `calc_se`); NULL (absent) otherwise,
  meaning analytic SEs are available and
  [`vcov`](https://rdrr.io/r/stats/vcov.html)/[`confint`](https://rdrr.io/r/stats/confint.html)/[`tidy`](https://netify-dev.github.io/sir/reference/tidy.sir.md)
  report the cluster-robust sandwich by default.

- dynamic_W:

  Logical, TRUE if W was time-varying (4D).

- symmetric/gamma/operator/rho_A/gain/stationary:

  Present for symmetric (A = B) fits: `symmetric = TRUE`,
  `operator = "symmetric"`; `gamma` is the full shared-influence vector
  (length p, all estimated, identified up to global sign,
  largest-magnitude gamma fixed positive); `rho_A` is the spectral
  radius of A (the max over periods for dynamic W),
  `gain = rho_A^2 / (n - 1)` is the stationarity gain, and `stationary`
  is `FALSE` when `gain >= 1` (the operator is explosive and estimates
  may be degenerate).

## Details

The SIR model specifies the expected outcome for the directed edge from
node i to node j at time t as:

\$\$g(\mu\_{i,j,t}) = \theta^T z\_{i,j,t} + \sum\_{k,l} X\_{k,l,t}
A\_{i,k} B\_{j,l}\$\$

Where:

- \\g(\cdot)\\ is the family link function (identity for normal, log for
  poisson, logit for binomial). The exogenous and bilinear terms enter
  the *linear predictor*, not the mean directly, so for poisson and
  binomial they act on the log / logit scale.

- \\\mu\_{i,j,t}\\ is the expected value of the outcome Y_ijt

- \\\theta\\ is a q-dimensional vector of coefficients for exogenous
  covariates

- \\z\_{i,j,t}\\ is a q-dimensional vector of exogenous covariates

- \\X\_{k,l,t}\\ represents the network state (often lagged Y) that
  carries influence

- \\A\_{i,k}\\ represents how node i is influenced by the behavior of
  node k

- \\B\_{j,l}\\ represents how node j's reception is affected by node l's
  position

The bilinear term \\\sum\_{k,l} X\_{k,l,t} A\_{i,k} B\_{j,l}\\ captures
network influence and can be parameterized using influence covariates W
through:

- \\A = \sum\_{r=1}^{p} \alpha_r W_r\\ (sender effects, \\\alpha_1 = 1\\
  fixed)

- \\B = \sum\_{r=1}^{p} \beta_r W_r\\ (receiver effects)

This parameterization reduces the number of parameters from \\O(m^2)\\
to \\O(p)\\, where \\p \ll m\\.

## Estimation Methods

**Alternating GLM/IRLS updates (method = "ALS"):**

- Iteratively updates the sender and receiver influence coefficients
  with GLM/IRLS subproblems

- Generally more stable for high-dimensional problems

- Better for sparse networks or when p is large

- May converge to local optima

**Direct Optimization (optim):**

- Uses BFGS to optimize all parameters simultaneously

- Can be faster for small problems

- May provide better solutions when good starting values are available

- More prone to numerical issues in high dimensions

Both engines report the same public parameterization with \\\alpha_1 =
1\\. BFGS optimizes that reduced vector directly; the default
alternating engine updates working alpha/beta coordinates with GLM/IRLS
subproblems and then normalizes to \\\alpha_1 = 1\\ when possible. The
full bilinear log-likelihood is often flat along a ridge: the
alternating engine and BFGS typically reach the same log-likelihood (to
well within 1%) yet can return parameter vectors that differ by a few
tenths on weakly identified influence coefficients. Compare fits by
log-likelihood / fitted values rather than by raw coefficients, and use
`fix_receiver = TRUE` to remove the alpha/beta scaling ambiguity when a
simpler sender-side channel is substantively adequate.

## Distribution Families

**Poisson:** For count data (e.g., number of interactions)

- Link function: log

- Variance function: \\V(\mu) = \mu\\

- Use when: Y_ijt represents counts

**Normal:** For continuous data (e.g., trade volumes, distances)

- Link function: identity

- Variance function: \\V(\mu) = \sigma^2\\

- Use when: Y_ijt is continuous and approximately normal

**Binomial:** For binary data (e.g., presence/absence of ties)

- Link function: logit

- Variance function: \\V(\mu) = \mu(1 - \mu)\\

- Use when: Y_ijt is binary (0/1)

## References

Minhas, S. & Hoff, P. D. (2025). Social Influence Regression. Political
Analysis.

## Examples

``` r
# \donttest{
set.seed(123)
m <- 8; T_len <- 5; p <- 2
Y <- array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
X <- array(0, dim = c(m, m, T_len))
X[,,2:T_len] <- Y[,,1:(T_len - 1)]
W <- array(rnorm(m * m * p), dim = c(m, m, p))
model <- sir(Y = Y, W = W, X = X, family = "poisson",
             method = "ALS", calc_se = FALSE, max_iter = 10)
print(model)
#> 
#> Social Influence Regression Model
#> 8 nodes, 5 time periods (directed)
#> Config: poisson | ALS
#> Status: converged | N = 280 | Log-Lik: -575.48 | AIC: 1157
#> Coefficients:
#>             Estimate
#> (alphaW) W2  -0.0912
#> (betaW) W1   -0.0037
#> (betaW) W2   -0.0066
#> (SEs not computed)
#> Use `summary()` for detailed results
coef(model)
#>  (alphaW) W2   (betaW) W1   (betaW) W2 
#> -0.091192198 -0.003650105 -0.006634265 
# }
```
