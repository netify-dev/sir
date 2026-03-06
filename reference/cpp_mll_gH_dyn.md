# Calculate Gradient and Hessian with Dynamic W

Computes gradient and Hessian of the negative log-likelihood when
influence covariates W vary over time (4D: m x m x p x T). The W array
is passed as a list of cubes to avoid Rcpp 4D array limitations.

## Usage

``` r
cpp_mll_gH_dyn(tab, Y, W_field, X, Z_list, family)
```

## Arguments

- tab:

  Parameter vector \[theta, alpha_2:p, beta\].

- Y:

  Three-dimensional array (m x m x T) of outcomes.

- W_field:

  List of T cubes, each m x m x p.

- X:

  Three-dimensional array (m x m x T).

- Z_list:

  List of q cubes (m x m x T), one per covariate.

- family:

  Distribution family string.

## Value

List with grad, hess, shess (after identifiability projection).
