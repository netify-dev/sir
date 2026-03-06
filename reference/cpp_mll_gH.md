# Calculate Gradient and Hessian for Direct Optimization

Computes the gradient vector and Hessian matrix of the negative
log-likelihood for the SIR model. These are essential for gradient-based
optimization methods like BFGS used in the direct optimization approach.

## Usage

``` r
cpp_mll_gH(tab, Y, W, X, Z_list, family)
```

## Arguments

- tab:

  Parameter vector \[theta, alpha_2:p, beta\] of length q+2p-1.

- Y:

  Three-dimensional array (m x m x T) of observed outcomes. Can contain
  NA values which are automatically skipped.

- W:

  Three-dimensional array (m x m x p) of influence covariates.

- X:

  Three-dimensional array (m x m x T) carrying network influence.

- Z_list:

  List of q three-dimensional arrays (m x m x T), one per covariate.
  Passed as list for efficient memory handling of 4D structure.

- family:

  String specifying distribution: "poisson", "normal", or "binomial".
  Determines the link function and variance structure.

## Value

List containing:

- grad:

  Gradient vector of length q+2p-1 (after identifiability projection)

- hess:

  Hessian matrix of dimension (q+2p-1) x (q+2p-1)

- shess:

  Score outer product matrix (for robust standard errors)

## Details

This function implements the analytical derivatives of the SIR
likelihood with respect to all parameters. The computation leverages the
bilinear structure of the model for efficiency.

**Gradient Computation:** For each parameter, the gradient accumulates:
\$\$\nabla\_{\cdot} NLL = \sum\_{i,j,t} (\mu\_{ijt} - y\_{ijt})
\frac{\partial \eta\_{ijt}}{\partial \cdot}\$\$

Where the partial derivatives are: - \\\partial \eta / \partial \theta_k
= Z\_{ijk,t}\\ - \\\partial \eta / \partial \alpha_k = \[W_k X
B'\]\_{ij}\\ - \\\partial \eta / \partial \beta_k = \[A X W_k'\]\_{ij}\\

**Hessian Computation:** The Hessian has two components: 1. Fisher
Information (always positive semi-definite): \$\$H\_{Fisher} =
\sum\_{i,j,t} w\_{ijt} \nabla \eta\_{ijt} \nabla \eta\_{ijt}'\$\$ where
\\w\_{ijt}\\ is the GLM weight (variance function)

2\. Observed Information adjustment (for non-canonical links):
\$\$H\_{Obs} = \sum\_{i,j,t} (\mu\_{ijt} - y\_{ijt}) \nabla^2
\eta\_{ijt}\$\$ This term captures the curvature from the bilinear
structure

**Identifiability Constraint:** Since alpha_1 is fixed at 1 for
identifiability, the function: - Computes full derivatives internally -
Projects out the alpha_1 dimension using matrix J - Returns
reduced-dimension gradient and Hessian

**Numerical Stability:** - Handles missing data (NA) by skipping those
observations - Adds small stabilization to binomial weights to prevent
singularity - Uses log-space computations where possible

**Computational Efficiency:** - Pre-computes A and B matrices once per
function call - Reuses matrix products across gradient components -
Vectorized operations using Armadillo - Avoids redundant calculations in
inner loops

## Note

The Hessian may not be positive definite far from the optimum, which can
cause optimization issues. The score outer product (shess) provides an
alternative for standard error calculation.

## Examples
