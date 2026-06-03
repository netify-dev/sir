# Get Scenario Values for Prediction

Computes representative covariate values for model-implied scenario
analysis. For each influence covariate, returns a set of values (the
10th, 25th, 50th, 75th, and 90th percentiles) at which to evaluate
predicted influence.

## Usage

``` r
get_scen_vals(
  data,
  vars = NULL,
  time = NULL,
  directed = TRUE,
  include_diag = FALSE,
  node_names = NULL,
  scen_vals = NULL
)
```

## Arguments

- data:

  A 3D (`m x m x p`) or 4D (`m x m x p x T`) array of influence
  covariates, or a data frame of covariates.

- vars:

  Character vector of variable names to compute scenario values for. If
  NULL, uses all available variables.

- time:

  Optional time index for time-varying covariates. If NULL, values are
  pooled across time.

- directed:

  Logical. If TRUE (default), treats the network as directed when
  computing scenario values. For square arrays, diagonals are excluded
  unless `include_diag = TRUE`; if FALSE, only the upper triangle is
  used.

- include_diag:

  Logical. If TRUE, include square-array diagonal cells when computing
  quantiles. Keep the default FALSE for one-mode networks; use TRUE for
  square bipartite sender/receiver arrays whose diagonal cells are real.

- node_names:

  Optional character vector of node names for labeling.

- scen_vals:

  Optional named list of pre-computed scenario values.

## Value

A named list mapping each variable to a vector of quantile values.

## Examples

``` r
dat <- sim_sir(m = 8, T_len = 10, p = 2, q = 1, family = "poisson", seed = 1)
sv <- get_scen_vals(dat$W)
str(sv)
#> List of 2
#>  $ var1: Named num [1:5] -1.202 -0.803 -0.198 0.453 1.039
#>   ..- attr(*, "names")= chr [1:5] "10%" "25%" "50%" "75%" ...
#>  $ var2: Named num [1:5] -1.043 -0.624 0.036 0.691 1.378
#>   ..- attr(*, "names")= chr [1:5] "10%" "25%" "50%" "75%" ...
```
