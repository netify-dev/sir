# Extract Fitted Values from a SIR Model

Returns the fitted values on the response scale: expected counts for
Poisson, probabilities for binomial, or conditional means for normal.

## Usage

``` r
# S3 method for class 'sir'
fitted(object, ...)
```

## Arguments

- object:

  A fitted `sir` object from
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md).

- ...:

  Additional arguments (unused).

## Value

An array with the same dimensions as `Y` (n1 x n2 x T) containing fitted
values on the response scale.
