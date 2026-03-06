# Print a Fitted SIR Model

Displays a compact overview of the fitted model: network dimensions,
family, method, convergence status, log-likelihood, AIC, and coefficient
estimates. Use [`summary()`](https://rdrr.io/r/base/summary.html) for a
more detailed report with p-values and influence matrix summaries.

## Usage

``` r
# S3 method for class 'sir'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  A fitted `sir` object from
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md).

- digits:

  Number of digits to print. Default uses `getOption("digits") - 3`.

- ...:

  Additional arguments (unused).

## Value

Invisibly returns the sir object.
