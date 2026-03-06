# Prepare Z Array for C++ Consumption

Converts a 4D array of exogenous covariates into a format optimized for
C++ processing. RcppArmadillo handles 4D arrays most efficiently when
decomposed into a list of 3D cubes.

## Usage

``` r
prepare_Z_list(Z)
```

## Arguments

- Z:

  Array of exogenous covariates. Can be:

  - NULL: Returns empty list

  - 3D array (m x m x T): Single time-varying covariate

  - 4D array (m x m x q x T): Multiple covariates

## Value

List of 3D arrays, where each element corresponds to one covariate
across all time periods. Length equals q (number of covariates).

## Details

This function handles the dimensional restructuring needed for
computation in the C++ backend. The transformation preserves the
covariate structure while enabling vectorized operations in Armadillo.

Input dimensions:

- 3D array (m x m x T): Interpreted as single covariate (q=1)

- 4D array (m x m x q x T): Multiple covariates

Output structure:

- List of length q

- Each element is an (m x m x T) array for one covariate

## Examples

``` r
if (FALSE) { # \dontrun{
# Single covariate
Z_single <- array(rnorm(10*10*5), dim=c(10,10,5))
Z_list <- prepare_Z_list(Z_single)
length(Z_list)  # Returns 1

# Multiple covariates  
Z_multi <- array(rnorm(10*10*3*5), dim=c(10,10,3,5))
Z_list <- prepare_Z_list(Z_multi)
length(Z_list)  # Returns 3
} # }
```
