# Inference and model comparison

Inference for bilinear models requires more care than for standard GLMs.
The bilinear structure in the SIR model can produce ill-conditioned
Hessians, which in turn yield unreliable Wald-based standard errors.
This vignette covers the inference tools available in the package:
classical and robust variance-covariance estimation, fixed-receiver
models that circumvent the identification problem, bootstrap standard
errors, and confidence intervals. For basic model fitting see
[`vignette("sir_overview")`](https://netify-dev.github.io/sir/articles/sir_overview.md).

## Setup

We use the same simulated data as in the overview vignette.

``` r
set.seed(42)
m = 15; T_len = 10; p = 2

W = array(0, dim = c(m, m, p))
geo = matrix(runif(m * m), m, m); geo = (geo + t(geo)) / 2; diag(geo) = 0
W[,,1] = geo
groups = sample(1:3, m, replace = TRUE)
W[,,2] = outer(groups, groups, "==") * 1.0; diag(W[,,2]) = 0
dimnames(W) = list(paste0("n", 1:m), paste0("n", 1:m), c("proximity", "shared_group"))

Z = array(0, dim = c(m, m, 1, T_len))
dist_mat = matrix(rnorm(m * m), m, m); dist_mat = (dist_mat + t(dist_mat)) / 2
diag(dist_mat) = NA
for (t in 1:T_len) Z[,,1,t] = dist_mat
dimnames(Z)[[3]] = "distance"

alpha_true = c(1, 0.3); beta_true = c(0.5, 0.2); theta_true = -0.1
A_true = alpha_true[1] * W[,,1] + alpha_true[2] * W[,,2]
B_true = beta_true[1] * W[,,1] + beta_true[2] * W[,,2]

Y = array(NA, dim = c(m, m, T_len))
Y[,,1] = matrix(rpois(m * m, 2), m, m); diag(Y[,,1]) = NA
for (t in 2:T_len) {
    X_t = log(Y[,,t-1] + 1); X_t[is.na(X_t)] = 0
    eta = theta_true * Z[,,1,t] + A_true %*% X_t %*% t(B_true)
    lambda = pmin(exp(eta), 100)
    Y[,,t] = matrix(rpois(m * m, c(lambda)), m, m); diag(Y[,,t]) = NA
}

X = array(0, dim = c(m, m, T_len))
for (t in 2:T_len) X[,,t] = log(Y[,,t-1] + 1)
X[is.na(X)] = 0
```

``` r
fit_als = sir(
    Y = Y, W = W, X = X, Z = Z,
    family = "poisson", method = "ALS", calc_se = TRUE
)
```

## Variance-covariance estimation

The package provides both classical (Hessian-based) and robust
(sandwich) variance-covariance estimates. The classical estimate
$H^{- 1}$ is efficient when the model is correctly specified. The
sandwich estimator $H^{- 1}SH^{- 1}$, where $S$ is the empirical score
covariance, remains valid under misspecification. In practice, the two
estimates often diverge for the bilinear parameters more than for the
exogenous covariate coefficients, which can itself serve as a diagnostic
for model adequacy.

``` r
V_classical = vcov(fit_als, type = "classical")
V_robust = vcov(fit_als, type = "robust")

# classical and robust standard errors
if (!is.null(V_classical)) sqrt(diag(V_classical))
#>          (Z) distance (alphaW) shared_group     (betaW) proximity 
#>          3.808729e-03          9.218215e-03          6.868159e-05 
#>  (betaW) shared_group 
#>          7.401193e-05
if (!is.null(V_robust)) sqrt(diag(V_robust))
#>          (Z) distance (alphaW) shared_group     (betaW) proximity 
#>          0.0280850001          0.0627795083          0.0004841839 
#>  (betaW) shared_group 
#>          0.0005015658
```

## Fixed-receiver model

Setting `fix_receiver = TRUE` constrains $\mathbf{B} = \mathbf{I}$ and
estimates only sender-side influence. This eliminates the bilinear
identification problem entirely and reduces the model to a standard GLM
with well-conditioned standard errors. When the research question
concerns sender-side influence alone, the fixed-receiver specification
is a natural starting point. It also provides a useful baseline for
comparison with the full bilinear model.

``` r
fit_fix = sir(
    Y = Y, W = W, X = X, Z = Z,
    family = "poisson",
    method = "ALS",
    fix_receiver = TRUE,
    calc_se = TRUE
)

summary(fit_fix)
#>                        Estimate Std. Error z value Pr(>|z|)    
#> (Z) distance          0.0121115  0.0034711   3.489 0.000484 ***
#> (alphaW) proximity    0.1159768  0.0003548 326.907  < 2e-16 ***
#> (alphaW) shared_group 0.0487102  0.0006104  79.795  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Model comparison

Information criteria provide a straightforward basis for comparing
nested specifications. The full bilinear model introduces additional
receiver-side parameters; whether the improvement in fit justifies the
added complexity is an empirical question.

``` r
AIC(fit_als)
#> [1] 201745.3
AIC(fit_fix)
#> [1] 170859.2
BIC(fit_als)
#> [1] 201767.9
BIC(fit_fix)
#> [1] 170876.1
```

## Bootstrap inference

When the Hessian is ill-conditioned, Wald-based standard errors become
unreliable. Bootstrap standard errors address this by constructing an
empirical sampling distribution for the parameter estimates. The block
bootstrap resamples entire time periods, preserving the temporal
dependence structure within each period. A parametric bootstrap, which
simulates new outcome arrays from the fitted model, is also available.

``` r
boot_results = boot_sir(
    fit_als,
    R = 100,
    type = "block",
    seed = 123,
    trace = FALSE
)

boot_results
#>                       Estimate Boot SE  2.5 % 97.5 %
#> (Z) distance            0.0465  0.0086 0.0270 0.0582
#> (alphaW) shared_group   0.5757  0.0107 0.5593 0.5982
#> (betaW) proximity       0.0109  0.0001 0.0109 0.0111
#> (betaW) shared_group    0.0062  0.0001 0.0060 0.0064
```

## Confidence intervals

Wald-based confidence intervals use the Hessian standard errors and
assume approximate normality of the sampling distribution.

``` r
confint(fit_als)
#>            2.5 %      97.5 %
#> [1,] 0.039065782 0.053995725
#> [2,] 0.557612194 0.593746934
#> [3,] 0.010814720 0.011083946
#> [4,] 0.006058405 0.006348527
```

When bootstrap results are available, percentile intervals constructed
from the empirical distribution of bootstrap replicates can be used
instead. These are generally more reliable for bilinear parameters,
particularly when the Hessian is ill-conditioned.

``` r
confint(fit_als, boot = boot_results)
#>           2.5 %      97.5 %
#> [1,] 0.02699795 0.058233021
#> [2,] 0.55934364 0.598212656
#> [3,] 0.01085486 0.011055276
#> [4,] 0.00600889 0.006419507
```
