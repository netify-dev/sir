# Flatten Y Array for GLM Input

Converts a 3D network outcome array into a vector suitable for GLM
estimation. Uses column-major ordering (R's default) to ensure
consistency with design matrices.

## Usage

``` r
flatten_Y(Y)
```

## Arguments

- Y:

  Three-dimensional array (m x m x T) of network outcomes.

## Value

Numeric vector of length m\*m\*T containing flattened outcomes.

## Details

The flattening follows R's column-major convention:

- Elements are ordered: Y\[1,1,1\], Y\[2,1,1\], ..., Y\[m,1,1\],
  Y\[1,2,1\], ...

- This matches how design matrices are constructed

- Preserves the correspondence between outcomes and covariates
