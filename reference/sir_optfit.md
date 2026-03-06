# Fit SIR Model via Direct Optimization

Optimizes all SIR model parameters simultaneously using BFGS with
analytical gradients computed via C++. An alternative to the ALS method
that can be faster for small problems or when good starting values are
available.

## Usage

``` r
sir_optfit(Y, W, X, Z, family, trace = 0, start = NULL)
```

## Arguments

- Y:

  Three-dimensional array (m x m x T) of network outcomes. Missing
  values (NA) are handled by excluding them from the likelihood. Large
  networks (m \> 100) may cause memory issues.

- W:

  Three-dimensional array (m x m x p) of influence covariates. Each
  slice parameterizes the influence matrices. If NULL, no network
  influence structure is included (reduces to standard GLM).

- X:

  Three-dimensional array (m x m x T) representing the network state
  carrying influence. Typically lagged Y. Required if W is provided.

- Z:

  Four-dimensional array (m x m x q x T) of exogenous covariates. These
  enter the model linearly without network interactions.

- family:

  Character string: "poisson", "normal", or "binomial". Determines the
  likelihood and link function used.

- trace:

  Integer controlling optimizer output:

  - 0: No output (default)

  - 1: Final convergence report

  - 2: Progress at each iteration

  - 3+: Detailed debugging information

- start:

  Optional numeric vector of starting values \[theta, alpha, beta\].
  Length must equal q + 2p. If NULL, uses smart initialization. Good
  starting values dramatically improve convergence.

## Value

A list with class "sir_optim_fit" containing:

- tab:

  Vector of optimized parameters \[theta, alpha, beta\]

- A:

  The m x m sender effects matrix

- B:

  The m x m receiver effects matrix

- convergence:

  Convergence code from optim (0 = success)

- message:

  Convergence message from optimizer

- iterations:

  Number of function evaluations

- value:

  Final negative log-likelihood

- hessian:

  Approximate Hessian at optimum (if requested)

- gradient:

  Final gradient (should be near zero)

## Details

This function treats the SIR negative log-likelihood as a single
optimization problem over all parameters (theta, alpha, beta). It uses
the BFGS quasi-Newton method with analytical gradients from the C++
backend.

If no starting values are provided, theta is initialized by fitting a
GLM ignoring network effects, and alpha/beta are set to small random
values.

Convergence code 0 indicates success; code 1 means the maximum number of
iterations was reached. Only supports static (3D) W.

## Examples

``` r
if (FALSE) { # \dontrun{
# Small network example
m <- 10
T <- 15
p <- 2

Y <- array(rpois(m*m*T, 2), dim=c(m,m,T))
X <- array(0, dim=c(m,m,T))
X[,,2:T] <- Y[,,1:(T-1)]
W <- array(rnorm(m*m*p, sd=0.1), dim=c(m,m,p))

# Fit with direct optimization
fit_optim <- sir_optfit(Y, W, X, Z=NULL,
                        family="poisson", trace=1)
                        
# Check convergence
if(fit_optim$convergence == 0) {
  cat("Optimization successful\n")
  print(fit_optim$tab)
}

# Compare with ALS
fit_als <- sir_alsfit(Y, W, X, Z=NULL,
                      family="poisson")

# Often similar results but different paths
} # }
```
