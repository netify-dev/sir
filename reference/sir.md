# Social Influence Regression (SIR) Model

Fits a Social Influence Regression model for network data with social
influence effects. The SIR model captures how network connections
influence outcomes through bilinear interaction terms, allowing for both
sender and receiver effects in directed networks.

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
  kron_mode = FALSE,
  ...
)
```

## Arguments

- Y:

  A three-dimensional array of dimensions (m x m x T) containing the
  network outcomes. Y\[i,j,t\] represents the directed outcome from node
  i to node j at time t. Can contain NA values for missing observations.
  The diagonal (self-loops) can be included or excluded depending on the
  application.

- W:

  Optional influence covariate array, either:

  - **3D array** (m x m x p): Static influence covariates. W\[i,j,r\]
    represents the r-th covariate for the edge from i to j. The same W
    is used for all time periods.

  - **4D array** (m x m x p x T): Dynamic (time-varying) influence
    covariates. W\[i,j,r,t\] allows the influence structure to change
    over time. Parameters (alpha, beta) are still estimated jointly
    across all periods, but the influence matrices A_t and B_t vary
    with t. Only ALS method is supported for 4D W.

  Common choices include graph Laplacians, geographic distance matrices,
  or node-level covariates expanded to edge-level. If NULL or p=0, the
  model uses only identity matrices (no network influence structure).

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

  Character string specifying the estimation method. Either "ALS"
  (Alternating Least Squares) or "optim" (direct optimization via BFGS).
  Default is "ALS" which is generally more stable.

- calc_se:

  Logical indicating whether to calculate standard errors for the
  parameters. Standard errors are computed using the observed
  information matrix. Setting to FALSE speeds up computation when
  uncertainty quantification is not needed.

- fix_receiver:

  Logical. If TRUE, fixes B = I (identity matrix) and estimates only
  (theta, alpha). This eliminates the bilinear identification problem
  (scaling ambiguity between A and B) by removing the receiver influence
  channel. The model becomes a standard GLM, yielding proper standard
  errors. Appropriate when receiver effects are negligible. Default is
  FALSE.

- symmetric:

  Logical. If TRUE, treats the network as undirected (symmetric). The
  function symmetrizes Y by averaging upper and lower triangles, uses
  only upper-triangle observations for fitting, and sets
  `fix_receiver = TRUE` (since sender/receiver distinction is
  meaningless for undirected networks). Default is FALSE.

- bipartite:

  Logical or NULL. Indicates whether the network is bipartite (senders
  and receivers are distinct node sets). If NULL (the default),
  bipartite status is inferred from Y: non-square arrays (n1 != n2) are
  treated as bipartite. Set to TRUE explicitly for square arrays where
  senders and receivers are nonetheless distinct populations. Setting
  FALSE on a non-square Y raises an error. Bipartite networks require
  `fix_receiver = TRUE`.

- kron_mode:

  Logical. If TRUE, replaces separate (alpha, beta) with a single p x p
  coefficient matrix C, where C\[r,s\] is the weight on W_r X W_s'. This
  is a general fix for the bilinear identification problem. Not yet
  implemented. Default is FALSE.

- ...:

  Additional arguments passed to the fitting functions:

  - `trace`: Logical or integer controlling output verbosity.

  - `tol`: Convergence tolerance for ALS (default 1e-8).

  - `max_iter`: Maximum ALS iterations (default 100).

## Value

An object of class `"sir"` with the following components:

- summ:

  Data frame of parameter estimates with columns `coef`, `se` (classical
  SE), `rse` (robust/sandwich SE), `t_se` (z-statistic using classical
  SE), `t_rse` (z-statistic using robust SE). Row names identify each
  parameter.

- A:

  Sender influence matrix. For static W: n1 x n1 matrix. For dynamic
  (4D) W: n1 x n1 x T array. Off-diagonal entry A\[i,k\] measures how
  much node k's behavior (via X) shapes node i's outgoing ties. Diagonal
  is set to zero.

- B:

  Receiver influence matrix. Same dimensions as A. Off-diagonal entry
  B\[j,l\] measures how node l's position shapes node j's incoming ties.
  Identity when `fix_receiver = TRUE`. Diagonal is zeroed.

- tab:

  Numeric vector of all estimated parameters in order: \[theta_1, ...,
  theta_q, alpha_2, ..., alpha_p, beta_1, ..., beta_p\]. When
  `fix_receiver = TRUE`: \[theta_1, ..., theta_q, alpha_1, ...,
  alpha_p\].

- theta:

  Coefficients for exogenous covariates Z (length q).

- alpha:

  Full alpha vector including the fixed alpha_1 = 1 (length p). When
  `fix_receiver = TRUE`, all alpha are free.

- beta:

  Coefficients for receiver influence covariates (length p). Empty when
  `fix_receiver = TRUE`.

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

  Logical, TRUE if the network is bipartite (n1 != n2).

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

## Details

The SIR model specifies the expected outcome for the directed edge from
node i to node j at time t as:

\$\$\mu\_{i,j,t} = \theta^T z\_{i,j,t} + \sum\_{k,l} X\_{k,l,t} A\_{i,k}
B\_{j,l}\$\$

Where:

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

**Alternating Least Squares (ALS):**

- Iteratively optimizes A given B, then B given A

- Generally more stable for high-dimensional problems

- Better for sparse networks or when p is large

- May converge to local optima

**Direct Optimization (optim):**

- Uses BFGS to optimize all parameters simultaneously

- Can be faster for small problems

- May provide better solutions when good starting values are available

- More prone to numerical issues in high dimensions

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
#> Status: converged | N = 320 | Log-Lik: -576.96 | AIC: 1159.9
#> Coefficients:
#>             Estimate
#> (alphaW) W2   3.9283
#> (betaW) W1    0.0010
#> (betaW) W2    0.0020
#> (SEs not computed)
#> Use `summary()` for detailed results
coef(model)
#> [1] 3.928314350 0.001010846 0.001972906
# }
```
