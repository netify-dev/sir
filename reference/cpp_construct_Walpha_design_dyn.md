# Construct Design Matrix for Beta Updates with Dynamic W

Builds the design matrix for updating receiver effects when influence
covariates vary over time (4D W array: m x m x p x T).

## Usage

``` r
cpp_construct_Walpha_design_dyn(W_field, X, alpha)
```

## Arguments

- W_field:

  List of T cubes, each m x m x p (one per time period).

- X:

  Three-dimensional array (m x m x T) of network states.

- alpha:

  Vector (p x 1) of current sender parameters.

## Value

Matrix (m\*m\*T x p) design matrix for beta GLM step.
