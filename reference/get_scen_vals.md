# Compute Summary Statistics for Scenario Construction

Extracts summary statistics (mean, quantiles) from a 4D covariate array.
These values are used to set baseline and varying levels when building
counterfactual scenario arrays with
[`get_scen_array`](https://netify-dev.github.io/sir/reference/get_scen_array.md).

## Usage

``` r
get_scen_vals(data, vars = NULL, time = NULL, directed = NULL)
```

## Arguments

- data:

  Four-dimensional array (m x m x p x T) of covariates, or a
  three-dimensional array (m x m x T) treated as a single variable.

- vars:

  Character vector of variable names (matching the third dimension of
  `data`) to summarise. If NULL, all variables are used.

- time:

  Integer vector of time indices to include. If NULL, all time periods
  are used.

- directed:

  Logical vector of length `length(vars)` indicating whether each
  variable is directed. For undirected variables the lower triangle is
  excluded before computing statistics. Defaults to TRUE for all
  variables.

## Value

A named list. Each element corresponds to one variable and contains:

- mean:

  Scalar mean of off-diagonal entries.

- quantiles:

  Named numeric vector of quantiles at 5% increments from 0 to 1.

## Examples

``` r
if (FALSE) { # \dontrun{
# W is m x m x p (or m x m x p x T)
vals <- get_scen_vals(W)
vals[["proximity"]]$mean
vals[["proximity"]]$quantiles
} # }
```
