# Getting started with sir

## Overview

Standard regression approaches to relational data treat each dyad as an
independent observation, ignoring the higher-order dependencies
generated when actors respond to one another’s past behavior. The `sir`
package fits Social Influence Regression (SIR) models that account for
these dependencies by regressing current network outcomes on the lagged
network state, decomposing influence into sender and receiver channels
parameterized by observable covariates.

The model specifies the expected outcome for directed edge $(i,j)$ at
time $t$ as:

$$\mu_{i,j,t} = {\mathbf{θ}}^{\top}\mathbf{z}_{i,j,t} + \sum\limits_{k,\ell}x_{k,\ell,t}\, a_{i,k}\, b_{j,\ell}$$

where $\mathbf{z}_{i,j,t}$ are exogenous covariates with coefficients
$\mathbf{θ}$, $x_{k,\ell,t}$ is the lagged network state (typically
$\log\left( y_{k,\ell,t - 1} + 1 \right)$), $a_{i,k}$ captures sender
influence (how predictive node $k$’s past sending is of node $i$’s
current sending), and $b_{j,\ell}$ captures receiver influence (how past
targeting of node $\ell$ predicts current targeting of node $j$). The
influence matrices $\mathbf{A}$ and $\mathbf{B}$ are parameterized
through influence covariates $\mathbf{W}$:

$$\mathbf{A} = \sum\limits_{r = 1}^{p}\alpha_{r}\mathbf{W}_{r}\qquad\mathbf{B} = \sum\limits_{r = 1}^{p}\beta_{r}\mathbf{W}_{r}$$

with $\alpha_{1} = 1$ fixed for identifiability. See the
[Methodology](https://netify-dev.github.io/sir/articles/methodology.md)
article for the full mathematical framework.

## Data format

The package expects data as multidimensional arrays:

| Array | Dimensions                     | Description                                                       |
|:------|:-------------------------------|:------------------------------------------------------------------|
| `Y`   | $m \times m \times T$          | Network outcomes (counts, continuous, or binary)                  |
| `X`   | $m \times m \times T$          | Network state carrying influence (typically lagged `Y`)           |
| `W`   | $m \times m \times p$          | Influence covariates parameterizing $\mathbf{A}$ and $\mathbf{B}$ |
| `Z`   | $m \times m \times q \times T$ | Exogenous dyadic covariates (optional)                            |

The
[`cast_array()`](https://netify-dev.github.io/sir/reference/cast_array.md)
utility converts long-format edge lists to the required array format
(see
[`vignette("sir_extensions")`](https://netify-dev.github.io/sir/articles/sir_extensions.md)).

## Simulated example data

To illustrate the workflow, we simulate a small directed network with
Poisson counts, an exogenous distance covariate, and two influence
covariates (geographic proximity and shared group membership). The
data-generating process follows the SIR specification directly, so we
know the true parameter values and can evaluate the model’s ability to
recover them.

``` r
set.seed(42)

m = 15     # number of nodes
T_len = 8  # number of time periods
p = 2      # number of influence covariates

# influence covariates
W = array(0, dim = c(m, m, p))

# W1: geographic proximity (symmetric, decays with distance)
geo = matrix(runif(m * m), m, m)
geo = (geo + t(geo)) / 2
diag(geo) = 0
W[,,1] = geo

# W2: shared group membership (block structure)
groups = sample(1:3, m, replace = TRUE)
shared = outer(groups, groups, "==") * 1.0
diag(shared) = 0
W[,,2] = shared

dimnames(W) = list(
    paste0("n", 1:m), paste0("n", 1:m),
    c("proximity", "shared_group")
)

# true influence parameters
alpha_true = c(1, 0.3)
beta_true  = c(0.5, 0.2)

# build true influence matrices
A_true = matrix(0, m, m)
for (r in 1:p) A_true = A_true + alpha_true[r] * W[,,r]
B_true = matrix(0, m, m)
for (r in 1:p) B_true = B_true + beta_true[r] * W[,,r]

# exogenous covariate: dyadic distance (time-invariant)
Z = array(0, dim = c(m, m, 1, T_len))
distance = matrix(rnorm(m * m, mean = 0, sd = 1), m, m)
distance = (distance + t(distance)) / 2
diag(distance) = NA
for (t in 1:T_len) Z[,,1,t] = distance
dimnames(Z)[[3]] = "distance"

theta_true = -0.1

# simulate Y as Poisson counts
Y = array(NA, dim = c(m, m, T_len))
dimnames(Y) = list(paste0("n", 1:m), paste0("n", 1:m), paste0("t", 1:T_len))

# first period: baseline rates
Y[,,1] = matrix(rpois(m * m, lambda = 2), m, m)
diag(Y[,,1]) = NA

# subsequent periods: influenced by lagged network
for (t in 2:T_len) {
    X_t = log(Y[,,t-1] + 1)
    X_t[is.na(X_t)] = 0
    eta = theta_true * Z[,,1,t] + A_true %*% X_t %*% t(B_true)
    lambda = exp(eta)
    lambda[lambda > 100] = 100
    Y[,,t] = matrix(rpois(m * m, lambda = c(lambda)), m, m)
    diag(Y[,,t]) = NA
}

# construct X as log-transformed lagged Y
X = array(0, dim = c(m, m, T_len))
for (t in 2:T_len) {
    X[,,t] = log(Y[,,t-1] + 1)
}
X[is.na(X)] = 0
```

This produces a network of 15 nodes over 8 periods. Edge counts are
Poisson, driven by a distance covariate ($\theta = - 0.1$) and network
influence flowing through geographic proximity and shared group
membership.

## Fitting the model

The default estimation method is Alternating Least Squares (ALS). The
algorithm iterates between GLM sub-problems for the sender and receiver
parameters until deviance converges. Each sub-problem is a standard GLM,
which makes the procedure numerically stable and efficient even for
moderately large networks.

``` r
fit_als = sir(
    Y = Y, W = W, X = X, Z = Z,
    family = "poisson",
    method = "ALS",
    calc_se = TRUE,
    trace = FALSE
)
```

### Model summary

``` r
summary(fit_als)
#>                        Estimate Std. Error z value Pr(>|z|)    
#> (Z) distance          4.591e-02  4.377e-03   10.49   <2e-16 ***
#> (alphaW) shared_group 5.667e-01  1.048e-02   54.09   <2e-16 ***
#> (betaW) proximity     1.095e-02  7.905e-05  138.58   <2e-16 ***
#> (betaW) shared_group  6.338e-03  8.573e-05   73.93   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Extracted coefficients

``` r
coef(fit_als)
#> [1] 0.045912767 0.566685792 0.010954690 0.006338205

logLik(fit_als)
#> 'log Lik.' -89985.73 (df=4)
AIC(fit_als)
#> [1] 179979.5
BIC(fit_als)
#> [1] 180001.2
```

### Influence matrices

The estimated $\mathbf{A}$ and $\mathbf{B}$ matrices contain sender and
receiver influence weights for each pair of nodes. Off-diagonal entries
indicate how predictive one actor’s past behavior is of another’s
current behavior. Examining the summary statistics of these matrices
provides a first indication of whether the model has recovered
meaningful influence structure.

``` r
# off-diagonal entries of the estimated influence matrices
A_hat = fit_als$A
A_offdiag = A_hat[row(A_hat) != col(A_hat)]
summary(A_offdiag)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.09552 0.46427 0.60420 0.68275 0.89105 1.50523

B_hat = fit_als$B
B_offdiag = B_hat[row(B_hat) != col(B_hat)]
summary(B_offdiag)
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> 0.001046 0.005086 0.006619 0.007518 0.009840 0.016620
```

## Diagnostic plots

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
produces a panel of diagnostics. The four default panels display
heatmaps of the influence matrices and the distribution of their
off-diagonal entries, which together indicate whether influence is
concentrated among particular dyads or distributed more broadly across
the network.

``` r
plot(fit_als, which = 1:4, title = "SIR Model Diagnostics")
```

![](sir_overview_files/figure-html/plot-diagnostics-1.png)

### Convergence history

For ALS fits, the deviance trajectory across iterations provides a
useful check that the algorithm has converged and that the final
estimates are not sensitive to the stopping point.

``` r
plot(fit_als, which = 5, combine = FALSE)
```

![](sir_overview_files/figure-html/plot-convergence-1.png)

### Coefficient plot

``` r
plot(fit_als, which = 6, combine = FALSE)
```

![](sir_overview_files/figure-html/plot-coefs-1.png)

## Fitted values and residuals

``` r
mu_hat = fitted(fit_als)
dim(mu_hat)
#> [1] 15 15  8
range(mu_hat, na.rm = TRUE)
#> [1]   0.9028613 668.7211119

r_dev = residuals(fit_als, type = "deviance")
r_prs = residuals(fit_als, type = "pearson")
r_rsp = residuals(fit_als, type = "response")
summary(as.vector(r_dev))
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
#> -21.13333  -0.03461   3.37267   4.52575   7.80506  28.73862       120
```

## Prediction

The [`predict()`](https://rdrr.io/r/stats/predict.html) method returns
values on either the link or response scale. For in-sample prediction,
the fitted model’s own data are used by default. For out-of-sample
prediction, new data can be supplied as a list.

``` r
pred_response = predict(fit_als, type = "response")
pred_link = predict(fit_als, type = "link")

range(pred_response, na.rm = TRUE)
#> [1]   0.9028613 668.7211119
range(pred_link, na.rm = TRUE)
#> [1] -0.1021863  6.5053671
```

``` r
new_pred = predict(
    fit_als,
    newdata = list(Y = Y_new, W = W, X = X_new, Z = Z_new),
    type = "response"
)
```

## Next steps

The remaining vignettes cover inference tools (variance-covariance
estimation, bootstrap standard errors, and confidence intervals) in
[`vignette("sir_inference")`](https://netify-dev.github.io/sir/articles/sir_inference.md),
and additional features (distribution families, symmetric and bipartite
networks, dynamic influence covariates, and data preparation utilities)
in
[`vignette("sir_extensions")`](https://netify-dev.github.io/sir/articles/sir_extensions.md).
The full mathematical framework is presented in
[`vignette("methodology")`](https://netify-dev.github.io/sir/articles/methodology.md).
