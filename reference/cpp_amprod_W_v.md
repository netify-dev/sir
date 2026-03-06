# Array-Matrix Product for Influence Matrices

Computes a weighted sum of influence covariate matrices to construct the
parameterized influence matrices A or B in the SIR model.

## Usage

``` r
cpp_amprod_W_v(W, v)
```

## Arguments

- W:

  Three-dimensional array (m x m x p) of influence covariates. Each
  slice W\[,,k\] represents one way nodes can influence each other
  (e.g., geographic proximity, social distance, shared attributes).

- v:

  Vector (p x 1) of coefficients for the linear combination. These are
  the parameters being estimated (either alpha or beta).

## Value

Matrix (m x m) representing the weighted combination of influence
covariate matrices. This becomes either the A or B matrix in the model.

## Details

This function implements the parameterization: Result = sum(k=1 to p)
v\[k\] \* W\[,,k\]

In the SIR model context: - A = sum(k=1 to p) alpha\[k\] \* W\[,,k\]
(alpha\[1\] = 1 fixed) - B = sum(k=1 to p) beta\[k\] \* W\[,,k\]

The parameterization reduces the number of free parameters from O(m^2)
to O(p), where typically p \<\< m. This makes estimation feasible for
larger networks.

Computational strategy: - Skips zero coefficients to save computation -
Uses in-place addition to minimize memory allocation - Leverages
Armadillo's expression templates for efficiency

## Note

The function checks for dimension compatibility and will throw an error
if v has incorrect length. Zero coefficients are detected and skipped to
improve performance when the model is sparse.

## Examples
