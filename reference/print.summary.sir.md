# Print a SIR Model Summary

Displays the full model summary including network dimensions, family,
method, a coefficient table with significance stars (when standard
errors are available), model fit statistics, convergence status, and
influence matrix summaries.

## Usage

``` r
# S3 method for class 'summary.sir'
print(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
)
```

## Arguments

- x:

  A `summary.sir` object from
  [`summary.sir`](https://netify-dev.github.io/sir/reference/summary.sir.md).

- digits:

  Number of significant digits to print. Default uses
  `getOption("digits") - 3`.

- signif.stars:

  Logical, whether to show significance stars beside p-values. Default
  uses `getOption("show.signif.stars")`.

- ...:

  Additional arguments (unused).

## Value

Invisibly returns the summary object.
