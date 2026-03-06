# Construct Design Matrix for Alpha Updates with Dynamic W

Builds the design matrix for updating sender effects when influence
covariates vary over time (4D W array: m x m x p x T).

## Usage

``` r
cpp_construct_Wbeta_design_dyn(W_field, X, beta)
```

## Arguments

- W_field:

  List of T cubes, each m x m x p (one per time period).

- X:

  Three-dimensional array (m x m x T) of network states.

- beta:

  Vector (p x 1) of current receiver parameters.

## Value

Matrix (m\*m\*T x p) design matrix for alpha GLM step.
