# Extract Number of Observations from a SIR Model

Returns the number of non-missing dyad-time observations used in
fitting. For square networks, diagonal entries (self-loops) are excluded
from this count.

## Usage

``` r
# S3 method for class 'sir'
nobs(object, ...)
```

## Arguments

- object:

  A fitted `sir` object from
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md).

- ...:

  Additional arguments (unused).

## Value

Integer count of observations.
