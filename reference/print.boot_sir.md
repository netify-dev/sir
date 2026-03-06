# Print Bootstrap SIR Results

Displays a table of point estimates, bootstrap standard errors, and 95%
percentile confidence intervals.

## Usage

``` r
# S3 method for class 'boot_sir'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  A `boot_sir` object from
  [`boot_sir`](https://netify-dev.github.io/sir/reference/boot_sir.md).

- digits:

  Number of significant digits. Default uses `getOption("digits") - 3`.

- ...:

  Additional arguments (unused).

## Value

Invisibly returns the `boot_sir` object.
