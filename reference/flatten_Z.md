# Flatten Z Array for GLM Input

Converts a 4D covariate array into a design matrix compatible with
flattened Y. Ensures proper alignment between outcomes and covariates
for GLM estimation.

## Usage

``` r
flatten_Z(Z)
```

## Arguments

- Z:

  Covariate array. Can be:

  - NULL: Returns NULL

  - 3D array (m x m x T): Converted to (m\*m\*T) x 1 matrix

  - 4D array (m x m x q x T): Converted to (m\*m\*T) x q matrix

## Value

Design matrix with dimensions (m\*m\*T) x q, or NULL if Z is NULL.
Column names are preserved from dimension names or auto-generated.

## Details

The flattening process maintains the correspondence between each
Y\[i,j,t\] and its associated covariates Z\[i,j,:,t\]. The resulting
matrix has:

- Rows: One for each outcome observation (m\*m\*T total)

- Columns: One for each covariate (q total)

- Order: Matches the flattening of Y (column-major)

Special handling for 3D input (single covariate):

- Automatically detected and converted to column matrix

- Preserves any dimension names for interpretability
