# Fit SIR Model via Alternating Least Squares (ALS)

Fits the SIR model by alternating between optimizing sender effects
(alpha) with receiver effects fixed, and vice versa. Each sub-step is a
standard GLM, making this approach more stable than direct optimization
for high-dimensional problems.

## Usage

``` r
sir_alsfit(
  Y,
  W,
  X,
  Z,
  family,
  trace = FALSE,
  tol = 1e-08,
  max_iter = 100,
  fix_receiver = FALSE,
  kron_mode = FALSE,
  dynamic_W = FALSE
)
```

## Arguments

- Y:

  Three-dimensional array (m x m x T) of network outcomes. Missing
  values (NA) are automatically handled by excluding them from the
  likelihood. The algorithm uses complete case analysis within each GLM
  step.

- W:

  Three-dimensional array (m x m x p) of influence covariates used to
  parameterize the influence matrices A and B. Each slice W\[,,r\]
  represents one influence covariate. If NULL or p=0, only identity
  matrices are used (no influence).

- X:

  Three-dimensional array (m x m x T) representing the network state
  that carries influence, typically lagged outcomes. Must be provided if
  W is non-NULL. This determines which network patterns affect future
  outcomes.

- Z:

  Four-dimensional array (m x m x q x T) of exogenous covariates, or
  NULL. These are covariates that directly affect outcomes but don't
  interact with the network influence structure. Examples include dyadic
  attributes or time trends.

- family:

  Character string specifying the GLM family: "poisson", "normal", or
  "binomial". This determines the link function and variance structure
  used in each GLM step.

- trace:

  Logical or integer controlling verbosity:

  - FALSE/0: No output

  - TRUE/1: Progress bar and convergence message

  - 2: Detailed iteration information including deviance

- tol:

  Numeric convergence tolerance. The algorithm stops when the relative
  change in deviance is less than this value. Default is 1e-8. Smaller
  values give more accurate results but require more iterations.

- max_iter:

  Integer maximum number of ALS iterations. Default is 100. Each
  iteration consists of one A-step and one B-step. Increase for
  difficult problems or when starting far from the optimum.

- fix_receiver:

  Logical. If TRUE, fixes B = I (identity matrix) and estimates only
  (theta, alpha) via a single GLM step. This eliminates the bilinear
  identification problem by removing the receiver influence channel. The
  model reduces to a standard GLM with well-conditioned standard errors.
  Default is FALSE.

- kron_mode:

  Logical. If TRUE, replaces separate (alpha, beta) with a single p x p
  coefficient matrix C. Not yet implemented. Default is FALSE.

- dynamic_W:

  Logical. If TRUE, W is treated as a 4D array (m x m x p x T) with
  time-varying influence covariates. Design matrices are constructed in
  R rather than C++. Default is FALSE.

## Value

A list with class "sir_als_fit" containing:

- tab:

  Vector of all parameters \[theta, alpha, beta\] in order

- A:

  The m x m sender effects matrix

- B:

  The m x m receiver effects matrix

- deviance:

  Final deviance (-2 \* log-likelihood + constant)

- iterations:

  Number of iterations until convergence

- converged:

  Logical indicating successful convergence

- THETA:

  Matrix tracking theta parameters across iterations

- ALPHA:

  Matrix tracking alpha parameters across iterations

- BETA:

  Matrix tracking beta parameters across iterations

- DEV:

  Matrix tracking deviance across iterations

- glm_alpha:

  Final GLM object from the A-step

- glm_beta:

  Final GLM object from the B-step

## Details

The algorithm exploits the bilinear structure: when B is fixed, the
model is linear in alpha (and theta), and vice versa. Each iteration
solves two GLMs using `speedglm` (if available) or
[`stats::glm`](https://rdrr.io/r/stats/glm.html).

Convergence is declared when the relative change in deviance falls below
`tol`, or when `max_iter` is reached. Supports both static (3D) and
dynamic (4D) W arrays.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example with simulated network data
m <- 15
T <- 20
p <- 2

# Generate network with influence
Y <- array(rpois(m*m*T, 3), dim=c(m,m,T))
X <- array(0, dim=c(m,m,T))
X[,,2:T] <- Y[,,1:(T-1)]  # Lagged Y

# Influence covariates (e.g., distance-based)
W <- array(rnorm(m*m*p), dim=c(m,m,p))

# Fit using ALS
fit <- sir_alsfit(Y, W, X, Z=NULL, family="poisson", 
                  trace=TRUE, tol=1e-6, max_iter=50)

# Examine convergence
plot(fit$DEV[,2], type="l", ylab="Deviance", xlab="Iteration")

# Extract influence matrices
A_matrix <- fit$A
B_matrix <- fit$B
} # }
```
