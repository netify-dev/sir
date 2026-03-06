# Construct Relational Covariates from a Network Array

Builds standard relational covariates from a base network array: the
original (main) effect, the reciprocal (transpose) effect, and a
transitive closure effect. These are common exogenous covariates (Z) in
SIR models capturing higher-order network dependencies.

## Usage

``` r
rel_covar(arr, name, effects = c("main", "reciprocal", "transitive"))
```

## Arguments

- arr:

  A 3D array (m x m x T) of network data. Typically an outcome variable
  or a covariate that varies across dyads and time.

- name:

  Character string used to label the covariates in the output array's
  third dimension. The main effect is labeled `name`, the reciprocal
  `paste0(name, "_recip")`, and the transitive `paste0(name, "_trans")`.

- effects:

  Character vector specifying which relational covariates to include.
  Any subset of `c("main", "reciprocal", "transitive")`. Default is all
  three.

## Value

A 4D array (m x m x q x T) where q is the number of requested effects (1
to 3). Suitable for passing directly as the `Z` argument to
[`sir`](https://netify-dev.github.io/sir/reference/sir.md).

## Details

- Main:

  The original array, Z_ij = arr_ij. This captures the direct dyadic
  effect.

- Reciprocal:

  The transpose, Z_ij = arr_ji. Captures whether the reverse
  relationship matters.

- Transitive:

  A measure of shared connectivity: Z_ij = (S network. Captures triadic
  closure and transitivity effects.

## Examples

``` r
if (FALSE) { # \dontrun{
# Build relational covariates from trade data
Z_trade <- rel_covar(trade_array, "trade")
dim(Z_trade)  # m x m x 3 x T

# Use only main and reciprocal effects
Z_simple <- rel_covar(trade_array, "trade", effects = c("main", "reciprocal"))

# Pass to sir() as exogenous covariates
fit <- sir(Y, W, X, Z = Z_trade, family = "poisson")
} # }
```
