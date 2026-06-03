# Build Scenario Array for Prediction

Constructs an artificial influence covariate grid for model-implied
scenario prediction, where one variable is set to a scenario value while
the others are held at the mean of their supplied scenario values. The
resulting array can be passed to
[`predict.sir`](https://netify-dev.github.io/sir/reference/predict.sir.md)
via `newdata`. The diagonal is set to zero by default for one-mode
scenario arrays so self-influence channels are not reintroduced. For
empirical counterfactuals, prefer copying the observed `W` array and
modifying a theoretically defined subset of cells.

## Usage

``` r
get_scen_array(
  var_to_vary,
  scen_vals,
  node_names,
  var_names,
  n_time = NULL,
  value = NULL,
  zero_diag = TRUE
)
```

## Arguments

- var_to_vary:

  Character name of the variable to vary.

- scen_vals:

  Named list of scenario values (from
  [`get_scen_vals`](https://netify-dev.github.io/sir/reference/get_scen_vals.md)).

- node_names:

  Character vector of node names.

- var_names:

  Character vector of all variable names in the array.

- n_time:

  Optional integer number of time periods. If NULL (default), returns a
  static 3D `n x n x p` array. Set this to return a 4D time-varying
  `n x n x p x n_time` array.

- value:

  Numeric scenario value for `var_to_vary`. If NULL, uses the mean of
  that variable's scenario values.

- zero_diag:

  Logical. If TRUE (default), set square-array diagonals to zero. Use
  FALSE for square bipartite sender-side arrays whose diagonal cells are
  real observations.

## Value

A static 3D scenario array by default, or a 4D array when `n_time` is
supplied.

## Examples

``` r
dat <- sim_sir(m = 8, T_len = 10, p = 2, q = 1, family = "poisson", seed = 1)
sv <- get_scen_vals(dat$W)
arr <- get_scen_array(names(sv)[1], scen_vals = sv,
                      node_names = paste0("n", 1:8), var_names = names(sv))
dim(arr)
#> [1] 8 8 2
```
