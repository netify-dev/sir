# Extract Residuals from a SIR Model

Returns residuals of the specified type. Response residuals are raw (Y -
fitted). Pearson residuals are standardized by the variance function.
Deviance residuals are signed square roots of the individual deviance
contributions and are most useful for diagnostic plots.

## Usage

``` r
# S3 method for class 'sir'
residuals(object, type = c("deviance", "pearson", "response"), ...)
```

## Arguments

- object:

  A fitted `sir` object from
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md).

- type:

  Character string specifying residual type: `"deviance"` (default),
  `"pearson"`, or `"response"`.

- ...:

  Additional arguments (unused).

## Value

An array with the same dimensions as `Y` containing the requested
residuals. Contains NA where `Y` is missing.
