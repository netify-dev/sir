# Build Counterfactual Scenario Array

Constructs a 4D array of covariate values for counterfactual prediction,
varying one variable across its empirical range while holding others at
their mean.

## Usage

``` r
get_scen_array(var_to_vary, scen_vals, node_names, var_names)
```

## Arguments

- var_to_vary:

  Character string naming the variable to vary.

- scen_vals:

  Named list produced by
  [`get_scen_vals`](https://netify-dev.github.io/sir/reference/get_scen_vals.md).

- node_names:

  Character vector of node names (used for the first two array
  dimensions).

- var_names:

  Character vector of all variable names (third dimension of the
  output).

## Value

A four-dimensional array of dimension `m x m x p x S`, where `S` is the
number of unique quantile values for the varying variable. Each slice
along the fourth dimension holds one scenario. Off-diagonal entries for
the varying variable take its quantile value, while all other variables
are set to their mean. Diagonals are zero.

## Examples

``` r
if (FALSE) { # \dontrun{
vals <- get_scen_vals(W)
scen <- get_scen_array("proximity", vals,
                       node_names = rownames(W),
                       var_names  = dimnames(W)[[3]])
} # }
```
