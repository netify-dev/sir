# Construct Design Matrix for Alpha Updates in ALS

Builds the design matrix for updating sender effects (alpha parameters)
in the Alternating Least Squares algorithm, holding receiver effects
(beta) fixed.

## Usage

``` r
cpp_construct_Wbeta_design(W, X, beta)
```

## Arguments

- W:

  Three-dimensional array (m x m x p) of influence covariates. These
  parameterize how influence flows through the network.

- X:

  Three-dimensional array (m x m x T) of network states over time.
  Typically contains lagged outcomes that carry influence forward.

- beta:

  Vector (p x 1) of current receiver effect parameters. These are held
  fixed while updating alpha in this ALS step.

## Value

Matrix (m\*m\*T x p) that serves as the design matrix for GLM
estimation. Each column corresponds to one influence covariate, rows
match vectorized Y.

## Details

In the ALS algorithm, when updating alpha with beta fixed, the model
becomes linear in alpha. The design matrix for this GLM sub-problem has
columns corresponding to each influence covariate.

For covariate k and observation (i,j,t), the design matrix element is:
\[W\[,,k\] \* X\[,,t\] \* B'\]\[i,j\]

Where B = sum_l beta\[l\] \* W\[,,l\] is the current receiver effects
matrix.

The algorithm: 1. Compute B from current beta and W 2. For each
covariate k: - Calculate W\[,,k\] \* X \* B' for all time points -
Flatten to match the vectorized Y 3. Combine into design matrix

This C++ implementation is 10-100x faster than the equivalent R code
using loops or apply functions, making ALS feasible for larger networks.

## Note

This function is called once per ALS iteration. The resulting matrix can
be large (m^2 \* T x p), so memory usage should be considered for big
networks.

## Examples
