# sir: Social Influence Regression Models

## Overview

The `sir` package implements Social Influence Regression (SIR) models for analyzing network data with social influence effects. The package provides efficient methods for fitting models with bilinear effects using either Alternating Least Squares (ALS) or direct optimization approaches.

## Features

- Support for multiple outcome distributions: Poisson, Normal, and Binomial
- Two fitting methods:
  - **ALS (Alternating Least Squares)**: Block coordinate descent approach (recommended for stability)
  - **Direct Optimization**: BFGS optimization with analytical gradients
- Efficient C++ backend using RcppArmadillo for fast computation
- Standard error calculation with both classical and robust (sandwich) estimators
- Handles high-dimensional network data with temporal dynamics

## Installation

You can install the development version of `sir` from GitHub:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install the sir package
devtools::install_github("username/sir")
```

Or install from local directory:

```r
devtools::install(".")
```

## Model Specification

The SIR model implements the following specification:

```
μ_{i,j,t} = θ^T z_{i,j,t} + α^T X_{i,j,t} β
```

Where:
- `Y`: (m × m × T) array of network outcomes
- `W`: (m × m × p) array of influence covariates
- `X`: (m × m × T) array of bilinear covariates (typically lagged outcomes)
- `Z`: (m × m × q × T) array of exogenous covariates
- `θ`: Coefficients for exogenous effects
- `α`, `β`: Coefficients for bilinear influence effects

## Basic Usage

```r
library(sir)

# Fit a SIR model with Poisson outcomes using ALS
model <- sir(Y = Y_array, 
             W = W_array, 
             X = X_array, 
             Z = Z_array,
             family = "poisson",
             method = "ALS")

# View model summary
print(model)

# Access influence matrices
A_matrix <- model$A  # Sender influence matrix
B_matrix <- model$B  # Receiver influence matrix

# Get coefficient estimates with standard errors
summary_table <- model$summ
```

## Advanced Options

### Using Direct Optimization

```r
# Fit using BFGS optimization instead of ALS
model_optim <- sir(Y = Y_array, 
                   W = W_array, 
                   X = X_array, 
                   Z = Z_array,
                   family = "normal",
                   method = "optim",
                   trace = 1)  # Show optimization progress
```

### Control ALS Algorithm

```r
# Customize ALS convergence criteria
model_als <- sir(Y = Y_array, 
                 W = W_array, 
                 X = X_array, 
                 Z = Z_array,
                 family = "binomial",
                 method = "ALS",
                 tol = 1e-8,      # Convergence tolerance
                 max_iter = 100,  # Maximum iterations
                 trace = TRUE)    # Print iteration details
```

## Authors

- **Shahryar Minhas** - *Primary Author* - [shahryar.minhas@gmail.com](mailto:shahryar.minhas@gmail.com)
- **Peter Hoff** - *Author* - [peter.hoff@duke.edu](mailto:peter.hoff@duke.edu)

## License

This project is licensed under the GPL (>= 3) License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this package in your research, please cite:

```
Minhas, S. and Hoff, P. (2025). sir: Social Influence Regression Models. 
R package version 0.1.0.
```

## Contributing

Please feel free to submit issues, feature requests, or pull requests through GitHub.

## Acknowledgments

This package uses:
- [Rcpp](https://www.rcpp.org/) and [RcppArmadillo](http://dirk.eddelbuettel.com/code/rcpp.armadillo.html) for C++ integration
- [speedglm](https://cran.r-project.org/package=speedglm) for fast GLM fitting (optional)