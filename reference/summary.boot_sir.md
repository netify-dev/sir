# Summary of Bootstrap SIR Results

Displays detailed bootstrap results including the coefficient table,
significance indicators (whether the 95% CI excludes zero), and
bootstrap distribution summaries (mean, sd, median).

## Usage

``` r
# S3 method for class 'boot_sir'
summary(object, ...)
```

## Arguments

- object:

  A `boot_sir` object from
  [`boot_sir`](https://netify-dev.github.io/sir/reference/boot_sir.md).

- ...:

  Additional arguments (unused).

## Value

Invisibly returns the `boot_sir` object.
