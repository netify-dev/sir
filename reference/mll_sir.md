# Calculate Negative Log-Likelihood for SIR Model

Computes the negative log-likelihood for the SIR model under the
specified distributional family. Diagonal entries (self-ties) and NA
values are excluded. Used as the objective function for parameter
estimation.

## Usage

``` r
mll_sir(tab, Y, W, X, Z, family, fix_receiver = FALSE)
```

## Arguments

- tab:

  Numeric vector of parameters \[theta, alpha_2:p, beta_1:p\].

- Y:

  Three-dimensional array (m x m x T) of observed outcomes. Can contain
  NA values which are automatically excluded.

- W:

  Three-dimensional array (m x m x p) of influence covariates, or NULL
  for no network influence.

- X:

  Three-dimensional array (m x m x T) carrying network influence.

- Z:

  Array of exogenous covariates (3D or 4D), or NULL.

- family:

  Character string specifying the distribution:

  - "poisson": Count data with log link

  - "normal": Continuous data with identity link (assumes sigma=1)

  - "binomial": Binary data with logit link

- fix_receiver:

  Logical. If TRUE, B is fixed to identity and tab is parsed as \[theta,
  alpha_1:p\]. Default FALSE.

## Value

Numeric scalar giving the negative log-likelihood. Lower values indicate
better fit. Used for optimization.

## Note

The normal family assumes unit variance (sigma=1) for simplicity. The
actual variance is estimated separately if needed.

## Examples

``` r
if (FALSE) { # \dontrun{
# Poisson example
Y <- array(rpois(1000, 2), dim=c(10,10,10))
nll <- mll_sir(tab, Y, W, X, Z, "poisson")

# Binomial example with missing data
Y_binary <- array(rbinom(1000, 1, 0.3), dim=c(10,10,10))
Y_binary[1,1,1] <- NA  # Missing value
nll <- mll_sir(tab, Y_binary, W, X, Z, "binomial")
} # }
```
