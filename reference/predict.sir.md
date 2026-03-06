# Predictions from a Fitted SIR Model

Generates predictions from a fitted SIR model for the training data or
for new data. Predictions can be on the link scale (linear predictor) or
the response scale (expected counts, probabilities, or means).

## Usage

``` r
# S3 method for class 'sir'
predict(object, newdata = NULL, type = c("response", "link"), ...)
```

## Arguments

- object:

  A fitted `sir` object from
  [`sir`](https://netify-dev.github.io/sir/reference/sir.md).

- newdata:

  Optional named list with components `W` (3D or 4D array), `X` (3D
  array), and/or `Z` (3D or 4D array) for counterfactual prediction.
  Dimensions must match the original fit. Any component not supplied is
  taken from the original fit. If NULL (default), returns predictions
  for the training data. Note: unlike many R predict methods, `newdata`
  is a list of arrays, not a data frame.

- type:

  Character string: `"link"` for linear predictor (eta) or `"response"`
  for expected values on the original scale. Default is `"response"`.

- ...:

  Additional arguments (unused).

## Value

An array (n1 x n2 x T) of predicted values on the requested scale.

## Details

For scenario (counterfactual) analysis, supply modified arrays in
`newdata`. For example, to see how the network would change if a
covariate increased by one unit, pass the modified Z array while keeping
W and X from the original fit.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- sir(Y, W, X, Z = Z, family = "poisson")

# In-sample fitted values
pred <- predict(model)

# Scenario: what if Z increases by 1 unit?
Z_shift <- Z + 1
pred_scenario <- predict(model, newdata = list(Z = Z_shift))

# Compare mean predictions
mean(pred, na.rm = TRUE)
mean(pred_scenario, na.rm = TRUE)
} # }
```
