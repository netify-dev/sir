# Calculate Linear Predictor (eta) for SIR Model

Computes the linear predictor \\\eta\_{i,j,t} = \theta^T z\_{i,j,t} +
\sum\_{k,l} X\_{k,l,t} A\_{i,k} B\_{j,l}\\, combining exogenous
covariate effects with the bilinear network influence term. Handles both
static (3D) and dynamic (4D) W arrays.

## Usage

``` r
eta_tab(tab, W, X, Z, fix_receiver = FALSE)
```

## Arguments

- tab:

  Numeric vector of parameters. For the standard model, ordered as
  \[theta, alpha_2:p, beta_1:p\] (length q + 2p - 1). For fix_receiver
  mode, ordered as \[theta, alpha_1:p\] (length q + p).

- W:

  Three-dimensional array (m x m x p) of influence covariates. Each
  slice W\[,,r\] parameterizes the influence structure. If NULL or p=0,
  no network influence is included.

- X:

  Three-dimensional array (m x m x T) carrying network influence.
  Typically lagged outcomes. Required even if W is NULL (can be zeros).

- Z:

  Array of exogenous covariates. Can be:

  - NULL or q=0: No exogenous effects

  - 3D array (m x m x T): Single covariate

  - 4D array (m x m x q x T): Multiple covariates

- fix_receiver:

  Logical. If TRUE, B is fixed to the identity matrix and tab is parsed
  as \[theta, alpha_1:p\] with all alpha free. Default FALSE.

## Value

Three-dimensional array (m x m x T) of linear predictors. Each element
eta\[i,j,t\] is the linear predictor for outcome Y\[i,j,t\] before
applying the link function.

## Examples

``` r
if (FALSE) { # \dontrun{
# Setup
m <- 10; T <- 5; p <- 2; q <- 1
W <- array(rnorm(m*m*p), dim=c(m,m,p))
X <- array(rnorm(m*m*T), dim=c(m,m,T))
Z <- array(rnorm(m*m*q*T), dim=c(m,m,q,T))

# Parameter vector
tab <- c(0.5,      # theta (q=1)
         0.2,      # alpha_2 (alpha_1=1 fixed)
         0.3, 0.4) # beta_1, beta_2
         
# Compute linear predictor
eta <- eta_tab(tab, W, X, Z)
dim(eta)  # Returns c(10, 10, 5)
} # }
```
