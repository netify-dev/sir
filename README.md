# sir: Social Influence Regression Models

## Overview

The `sir` package implements Social Influence Regression models for network analysis where outcomes on directed edges are influenced by network structure itself. This framework addresses a challenge in network analysis: the interdependence between network outcomes and network structure. Regression approaches fail to account for the fact that observations in network data are dependent due to shared nodes and structural patterns.

## Motivation

Network data presents inferential challenges because the network structure that connects observations also creates dependencies between them. Consider international trade networks, social media interactions, or neural connectivity patterns - in each case, the outcome of interest (trade volume, message frequency, signal strength) both depends on and influences the network structure.

The Social Influence Regression (SIR) model provides an approach to this problem by modeling how network outcomes are influenced by the network itself through bilinear effects. The model decomposes network influence into sender effects (how nodes influence others) and receiver effects (how nodes are influenced by others), providing interpretable parameters while maintaining tractability.

## Model Specification

The SIR model specifies the expected outcome for directed edge $(i,j)$ at time $t$ as:

$$\mu_{i,j,t} = \mathbf{z}_{i,j,t}^T \boldsymbol{\theta} + \sum_{k,\ell} x_{k,\ell,t} a_{i,k} b_{j,\ell}$$

where:
- $\mu_{i,j,t}$ is the expected value of outcome $y_{i,j,t}$
- $\mathbf{z}_{i,j,t}$ is a vector of exogenous covariates with coefficients $\boldsymbol{\theta}$
- $x_{k,\ell,t}$ represents the network state carrying influence (often lagged outcomes)
- $a_{i,k}$ captures how node $i$ is influenced by node $k$'s sending behavior
- $b_{j,\ell}$ captures how node $j$'s receiving is affected by node $\ell$'s receiving behavior

The influence matrices $\mathbf{A}$ and $\mathbf{B}$ are parameterized as:

$$\mathbf{A} = \mathbf{I} + \sum_{r=1}^{p-1} \alpha_r \mathbf{W}_r$$
$$\mathbf{B} = \beta_0 \mathbf{I} + \sum_{r=1}^{p} \beta_r \mathbf{W}_r$$

where $\mathbf{W}_r$ are known influence covariates (e.g., geographic distance, shared attributes).

## Installation

### Install from GitHub Releases

Pre-built binaries are available. Download the appropriate file from the [releases page](https://github.com/YOUR_USERNAME/sir/releases) and install:

```r
# Windows
install.packages("sir_0.1.0_windows.zip", repos = NULL, type = "binary")

# macOS (Apple Silicon)
install.packages("sir_0.1.0_macos-arm64.tgz", repos = NULL, type = "binary")

# macOS (Intel)
install.packages("sir_0.1.0_macos-intel.tgz", repos = NULL, type = "binary")

# Linux (source)
install.packages("sir_0.1.0_linux.tar.gz", repos = NULL, type = "source")
```

### Install from GitHub (Development Version)

```r
# Install remotes if needed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install sir from GitHub
remotes::install_github("YOUR_USERNAME/sir")
```

### Install from Source

```r
# Clone the repository first, then:
install.packages(".", repos = NULL, type = "source")

# Or using devtools
devtools::install(".")
```


## Basic Usage

```r
library(sir)

# Prepare data as 3D/4D arrays
# Y: m × m × T array of network outcomes
# W: m × m × p array of influence covariates  
# X: m × m × T array (typically lagged Y)
# Z: m × m × q × T array of exogenous covariates

# Fit model using Alternating Least Squares
model <- sir(
  Y = Y,
  W = W,
  X = X,
  Z = Z,
  family = "poisson",
  method = "ALS",
  trace = TRUE
)

# Model summary with significance tests
summary(model)

# Diagnostic plots
plot(model)

# Extract components
coef(model)          # Coefficient estimates
model$A              # Sender effects matrix
model$B              # Receiver effects matrix
logLik(model)        # Log-likelihood
AIC(model)           # Information criteria
```

## Features

### Distribution Families

The package supports three distributional families:

```r
# Count data (e.g., number of interactions)
model_poisson <- sir(Y, W, X, Z, family = "poisson")

# Continuous data (e.g., trade volumes)
model_normal <- sir(Y, W, X, Z, family = "normal")

# Binary data (e.g., presence/absence of ties)
model_binomial <- sir(Y, W, X, Z, family = "binomial")
```

### Estimation Methods

Two estimation approaches are available:

```r
# Alternating Least Squares (recommended for stability)
model_als <- sir(Y, W, X, Z, family = "poisson", method = "ALS",
                 tol = 1e-8, max_iter = 100)

# Direct optimization via BFGS (can be faster for small problems)
model_optim <- sir(Y, W, X, Z, family = "poisson", method = "optim",
                   trace = 1)
```


### S3 Methods

The package provides S3 methods for model inspection:

```r
print(model)         # Basic model information
summary(model)       # Detailed summary with tests
coef(model)          # Extract coefficients
fitted(model)        # Fitted values
residuals(model)     # Residuals (various types)
logLik(model)        # Log-likelihood
AIC(model)           # Akaike Information Criterion
BIC(model)           # Bayesian Information Criterion
vcov(model)          # Variance-covariance matrix
predict(model, newdata)  # Predictions
plot(model)          # Diagnostic plots (ggplot2)
```

## Examples

### International Trade Network

```r
# Analyzing bilateral trade flows
# Y: trade volumes between countries over time
# W: inverse geographic distance, shared language indicators
# X: lagged trade volumes
# Z: GDP, population, trade agreements

trade_model <- sir(
  Y = trade_volumes,
  W = influence_matrices,
  X = lagged_trade,
  Z = country_covariates,
  family = "normal",
  method = "ALS"
)

# Identify influential trading partners
influential_exporters <- which(rowSums(abs(trade_model$A)) > threshold)
receptive_importers <- which(colSums(abs(trade_model$B)) > threshold)
```

## Mathematical Details

### Likelihood

The log-likelihood for the SIR model under the exponential family is:

$$\ell(\boldsymbol{\theta}, \boldsymbol{\alpha}, \boldsymbol{\beta}) = \sum_{i,j,t} \left[ y_{i,j,t} \eta_{i,j,t} - b(\eta_{i,j,t}) \right]$$

where $\eta_{i,j,t}$ is the linear predictor and $b(\cdot)$ is the log-partition function.

### Identifiability

For identifiability, we impose the constraint $\alpha_1 = 1$, fixing the scale of the bilinear term. This is automatically handled in the estimation procedures.

### Computational Complexity

- **ALS**: $O(Tm^2p + Tm^3)$ per iteration
- **Direct optimization**: $O(Tm^4)$ per gradient evaluation
- Memory requirement: $O(m^2(T + p))$

## Citation

If you use this package in your research, please cite:

```bibtex
@article{minhas:hoff:2025,
  title={Decomposing Network Dynamics: Social Influence Regression},
  author={Minhas, Shahryar and Hoff, Peter},
  journal={Political Analysis},
  year={2025},
}

@Manual{sir-package,
  title = {sir: Social Influence Regression Models},
  author = {Shahryar Minhas and Peter Hoff},
  year = {2025},
  note = {R package version 0.1.0},
  url = {https://github.com/YOUR_USERNAME/sir}
}
```

## Authors

- **Shahryar Minhas** (Michigan State University)
- **Peter Hoff** (Duke University)

## License

MIT License - see [LICENSE](LICENSE) file for details.
