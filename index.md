# sir: Social Influence Regression Models

## Motivation

Relational outcomes in political and social systems are rarely
independent. When country $i$ initiates conflict with country $k$ at
time $t$, that action may reshape the likelihood that country $j$ does
the same at time $t + 1$, particularly if $i$ and $j$ share an alliance,
geographic proximity, or a common adversary. Standard regression
approaches to network data treat each directed dyad as an independent
observation, conditioning on covariates but discarding the higher-order
dependencies that constitute network influence.

The Social Influence Regression (SIR) model addresses this limitation
directly. Rather than treating network influence as a nuisance or a
latent quantity to be absorbed by random effects, the model estimates
influence as a function of observable covariates. The framework operates
on longitudinal network data (a time series of $n \times n$ relational
matrices) and asks how past interactions across the network predict
current outcomes, and what dyad-level or actor-level covariates account
for those predictive patterns. The methodological framework is
introduced in Minhas and Hoff (2025), “Decomposing Network Dynamics:
Social Influence Regression,” *Political Analysis*.

## Model specification

The SIR model specifies the expected outcome for directed edge $(i,j)$
at time $t$ as:

$$\mu_{i,j,t} = \mathbf{z}_{i,j,t}^{T}{\mathbf{θ}} + \sum\limits_{k,\ell}x_{k,\ell,t}a_{i,k}b_{j,\ell}$$

The first term is a standard regression of the outcome on exogenous
covariates $\mathbf{z}_{i,j,t}$ with coefficients $\mathbf{θ}$. The
second term captures network influence: $a_{i,k}$ measures how
predictive node $k$’s past sending behavior is of node $i$’s current
sending, $b_{j,\ell}$ captures an analogous relationship on the receiver
side, and the network state $x_{k,\ell,t}$ (typically lagged outcomes)
carries the influence signal through these matrices.

The influence matrices $\mathbf{A}$ and $\mathbf{B}$ are parameterized
through known influence covariates $\mathbf{W}_{r}$ (geographic
distance, alliance ties, shared attributes):

$$\mathbf{A} = \sum\limits_{r = 1}^{p}\alpha_{r}\mathbf{W}_{r},\quad\alpha_{1} = 1{\mspace{6mu}\text{(fixed for identifiability)}}$$$$\mathbf{B} = \sum\limits_{r = 1}^{p}\beta_{r}\mathbf{W}_{r}$$

This parameterization reduces the parameter count from
$O\left( n^{2} \right)$ to $O(p)$, where $p$ is the number of influence
covariates. More importantly, it makes the estimated influence patterns
directly interpretable: the $\alpha$ and $\beta$ coefficients identify
which covariates matter for influence and by how much.

## Installation

### From GitHub releases

Pre-built binaries are available from the [releases
page](https://github.com/netify-dev/sir/releases):

``` r
# windows
install.packages("sir_0.1.0_windows.zip", repos = NULL, type = "binary")

# macOS (apple silicon)
install.packages("sir_0.1.0_macos-arm64.tgz", repos = NULL, type = "binary")

# macOS (intel)
install.packages("sir_0.1.0_macos-intel.tgz", repos = NULL, type = "binary")

# linux (source)
install.packages("sir_0.1.0_linux.tar.gz", repos = NULL, type = "source")
```

### From GitHub (development version)

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}
remotes::install_github("netify-dev/sir")
```

### From source

``` r
install.packages(".", repos = NULL, type = "source")
# or
devtools::install(".")
```

## Usage

``` r
library(sir)

# data inputs:
# Y: m x m x T array of network outcomes
# W: m x m x p array of influence covariates
# X: m x m x T array of lagged network state
# Z: m x m x q x T array of exogenous covariates

model = sir(
    Y = Y, W = W, X = X, Z = Z,
    family = "poisson",
    method = "ALS",
    trace = TRUE
)

summary(model)
plot(model)

# extract components
coef(model)          # parameter estimates
model$A              # sender influence matrix
model$B              # receiver influence matrix
confint(model)       # confidence intervals
```

## Estimation

The package provides two estimation approaches. The default, Alternating
Least Squares (ALS), exploits the bilinear structure of the model by
iterating between GLM sub-problems for $({\mathbf{θ}},{\mathbf{α}})$ and
$({\mathbf{θ}},{\mathbf{β}})$. Each sub-problem is a standard
generalized linear model, so the full estimation reduces to a sequence
of low-dimensional optimizations. This is generally more stable for
high-dimensional problems and substantially faster than the Bayesian
approach originally used for bilinear network autoregressions. A direct
BFGS optimizer (`method = "optim"`) is also available and uses
analytical gradients computed via C++; it can converge faster for small
networks but tends to be less stable when the number of influence
covariates is large.

Standard errors are computed from the Hessian (classical) and the
sandwich estimator (robust). When the Hessian is ill-conditioned, which
is common in bilinear models,
[`boot_sir()`](https://netify-dev.github.io/sir/reference/boot_sir.md)
provides bootstrap standard errors via block resampling of time periods
or parametric simulation.

## Distribution families

The model supports Poisson (count outcomes with log link), Normal
(continuous outcomes with identity link), and Binomial (binary outcomes
with logit link) families.

``` r
model_poisson  = sir(Y, W, X, Z, family = "poisson")
model_normal   = sir(Y, W, X, Z, family = "normal")
model_binomial = sir(Y, W, X, Z, family = "binomial")
```

## S3 methods

``` r
summary(model)                  # coefficient table with standard errors
coef(model)                     # parameter estimates
fitted(model)                   # fitted values on response scale
residuals(model, type = "deviance")  # deviance residuals
logLik(model)                   # log-likelihood
AIC(model); BIC(model)          # information criteria
vcov(model, type = "robust")    # sandwich variance-covariance
confint(model)                  # Wald confidence intervals
predict(model, newdata, type = "response")
plot(model, which = 1:6)        # diagnostic plots
```

## Additional features

The package supports symmetric (undirected) networks via
`symmetric = TRUE`, bipartite (rectangular) networks via
`bipartite = TRUE` or automatic detection from non-square $\mathbf{Y}$,
dynamic (time-varying) influence covariates via 4D `W` arrays,
counterfactual scenario construction via
[`get_scen_vals()`](https://netify-dev.github.io/sir/reference/get_scen_vals.md)
and
[`get_scen_array()`](https://netify-dev.github.io/sir/reference/get_scen_array.md),
and the
[`cast_array()`](https://netify-dev.github.io/sir/reference/cast_array.md)
utility for converting edge-list data to the required array format.

## Citation

``` bibtex
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
  url = {https://github.com/netify-dev/sir}
}
```

## Authors

Shahryar Minhas (Michigan State University) and Peter Hoff (Duke
University).

## License

MIT License. See [LICENSE](https://netify-dev.github.io/sir/LICENSE) for
details.
