# Construct Design Matrix for Beta Updates in ALS

Builds the design matrix for updating receiver effects (beta parameters)
in the Alternating Least Squares algorithm, holding sender effects
(alpha) fixed.

## Usage

``` r
cpp_construct_Walpha_design(W, X, alpha)
```

## Arguments

- W:

  Three-dimensional array (m x m x p) of influence covariates.

- X:

  Three-dimensional array (m x m x T) of network states over time.

- alpha:

  Vector (p x 1) of current sender effect parameters. These are held
  fixed while updating beta in this ALS step.

## Value

Matrix (m\*m\*T x p) serving as the design matrix for beta GLM
estimation.

## Details

In the ALS algorithm, when updating beta with alpha fixed, the model
becomes linear in beta. This function constructs the required design
matrix efficiently.

For covariate k and observation (i,j,t), the design matrix element is:
\[A \* X\[,,t\] \* W\[,,k\]'\]\[i,j\]

Where A = sum_l alpha\[l\] \* W\[,,l\] is the current sender effects
matrix (with alpha\[1\] = 1 fixed for identifiability).

The algorithm mirrors the alpha update but with roles reversed: 1.
Compute A from current alpha and W 2. For each covariate k: - Calculate
A \* X \* W\[,,k\]' for all time points - Flatten to match the
vectorized Y 3. Combine into design matrix

## Note

The symmetry with cpp_construct_Wbeta_design reflects the bilinear
structure of the model, where sender and receiver effects play dual
roles.

## Examples
