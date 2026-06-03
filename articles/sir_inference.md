# Inference and model comparison

Inference for bilinear models takes more care than for an ordinary GLM.
The bilinear influence term can produce an ill-conditioned Hessian,
which makes naive Wald standard errors unreliable. This vignette covers
the tools the package provides: classical and robust variance estimates,
confidence intervals, bootstrap and jackknife checks, and
information-criterion model comparison. For basic fitting see
[`vignette("sir_overview")`](https://netify-dev.github.io/sir/articles/sir_overview.md).

| Tool | What it is useful for | Main limitation |
|:---|:---|:---|
| Classical SE | Fast model-based Wald intervals | Assumes independent dyad-period scores and a stable Hessian |
| HC0 robust SE | Heteroskedasticity or overdispersion | Does not address shared-actor network dependence |
| Multiway cluster SE | Sender, receiver, and time dependence in the score | Still uses the Hessian as bread, so weak identification remains a problem |
| Block bootstrap | Time-period resampling | Preserves within-period dependence but not actor dependence |
| Dyad jackknife | Delete-one-actor sensitivity and actor dependence checks | Normal jackknife intervals, not bootstrap percentile intervals |

## Setup

We simulate a correctly specified Poisson model with
[`sim_sir()`](https://netify-dev.github.io/sir/reference/sim_sir.md) so
the standard errors are well behaved and we can focus on the mechanics.

``` r

set.seed(7)
dat <- sim_sir(m = 14, T_len = 80, p = 2, q = 2, family = "poisson", seed = 7)
fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z,
           family = "poisson", calc_se = TRUE, seed = 1)
```

## Classical (Hessian-based) standard errors

`sir(calc_se = TRUE)` computes classical standard errors from the
inverse observed-information matrix. These are what
[`summary()`](https://rdrr.io/r/base/summary.html) prints for quick
model inspection. For reporting, the default accessors now use
cluster-robust uncertainty when the fit supports it; request classical
intervals explicitly when you need a Hessian-only comparison:

``` r

sqrt(diag(vcov(fit, type = "classical")))   # the SEs reported by summary()
#>      (Z) Z1      (Z) Z2 (alphaW) W2  (betaW) W1  (betaW) W2 
#> 0.007030010 0.006552255 0.026505862 0.005033336 0.007607866
confint(fit, se.type = "classical")         # Hessian-only Wald intervals
#>                  2.5 %     97.5 %
#> (Z) Z1      -0.1807515 -0.1531944
#> (Z) Z2       0.4607029  0.4863873
#> (alphaW) W2  0.6826715  0.7865726
#> (betaW) W1  -0.3599133 -0.3401830
#> (betaW) W2  -0.2176376 -0.1878153
```

## Robust (sandwich) standard errors

The HC0 sandwich estimator $`H^{-1} S H^{-1}`$ (with $`S`$ the empirical
score covariance) is useful for heteroskedasticity or overdispersion
when dyad-period score contributions are otherwise independent and the
Hessian bread is stable. Compare it against the classical SEs:

``` r

se_classical <- sqrt(diag(vcov(fit, type = "classical")))
se_robust    <- sqrt(diag(vcov(fit, type = "robust")))

data.frame(
    term      = names(coef(fit)),
    classical = round(unname(se_classical), 4),
    robust    = round(unname(se_robust), 4),
    ratio     = round(unname(se_robust / se_classical), 2),
    row.names = NULL
)
#>          term classical robust ratio
#> 1      (Z) Z1    0.0070 0.0070  0.99
#> 2      (Z) Z2    0.0066 0.0065  0.99
#> 3 (alphaW) W2    0.0265 0.0261  0.98
#> 4  (betaW) W1    0.0050 0.0049  0.96
#> 5  (betaW) W2    0.0076 0.0075  0.99

confint(fit, se.type = "robust")
#>                  2.5 %     97.5 %
#> (Z) Z1      -0.1806005 -0.1533454
#> (Z) Z2       0.4608553  0.4862350
#> (alphaW) W2  0.6834737  0.7857704
#> (betaW) W1  -0.3595586 -0.3405378
#> (betaW) W2  -0.2174425 -0.1880104
```

For well-specified data the HC0 robust SEs barely move (ratio near 1):
HC0 corrects for overdispersion, not for the dyadic dependence that
dominates relational data. That dependence is what the
**cluster-robust** estimator addresses.

## Multiway cluster-robust standard errors

The $`(i, j)`$ and $`(j, i)`$ dyads share both actors, and every tie
involving an actor is correlated.
[`vcov()`](https://rdrr.io/r/stats/vcov.html) clusters the score
contributions on the sender, receiver, and time margins (`"twoway"`
remains an alias for this multiway estimator). For supported static
directed-network fits, report this interval in place of classical
Hessian SEs when the Hessian bread is stable:

``` r

se_cluster <- sqrt(diag(vcov(fit)))

data.frame(
    term      = names(coef(fit)),
    classical = round(unname(se_classical), 4),
    cluster   = round(unname(se_cluster), 4),
    ratio     = round(unname(se_cluster / se_classical), 2),
    row.names = NULL
)
#>          term classical cluster ratio
#> 1      (Z) Z1    0.0070  0.0062  0.88
#> 2      (Z) Z2    0.0066  0.0082  1.26
#> 3 (alphaW) W2    0.0265  0.0235  0.89
#> 4  (betaW) W1    0.0050  0.0047  0.93
#> 5  (betaW) W2    0.0076  0.0055  0.72
```

The cluster-robust SEs can move in either direction in this clean
simulation. In dependent real data they are often more conservative
because they stop assuming independent dyad-times. `confint(fit)`
returns the matching Wald intervals:

``` r

confint(fit)
#>                  2.5 %     97.5 %
#> (Z) Z1      -0.1790843 -0.1548617
#> (Z) Z2       0.4574220  0.4896683
#> (alphaW) W2  0.6885972  0.7806470
#> (betaW) W1  -0.3591747 -0.3409217
#> (betaW) W2  -0.2135366 -0.1919164
```

## Bootstrap inference

[`boot_sir()`](https://netify-dev.github.io/sir/reference/boot_sir.md)
refits the model on resampled or reduced data. The **block** bootstrap
resamples whole time periods (preserving within-period dependence); the
**dyad jackknife** deletes each actor in turn; the **parametric**
bootstrap simulates new outcome arrays from the fitted model.

``` r

br <- boot_sir(fit, R = 20, type = "block", seed = 123, trace = FALSE)
br_dyad <- boot_sir(fit, type = "dyad", seed = 123, trace = FALSE)

data.frame(
    method = c("Block Bootstrap", "Dyad Jackknife"),
    valid_refits = c(br$n_valid, br_dyad$n_valid),
    total_refits = c(br$n_total, br_dyad$n_total),
    interval_type = c(br$interval, br_dyad$interval)
)
#>            method valid_refits total_refits    interval_type
#> 1 Block Bootstrap           20           20       percentile
#> 2  Dyad Jackknife           14           14 normal-jackknife
```

Percentile intervals from block or parametric bootstrap replicates are
often more reliable than Wald intervals for the bilinear parameters:

``` r

confint(fit, boot = br)
#>                  2.5 %     97.5 %
#> (Z) Z1      -0.1822541 -0.1559235
#> (Z) Z2       0.4551020  0.4867091
#> (alphaW) W2  0.6766893  0.7749009
#> (betaW) W1  -0.3612847 -0.3429720
#> (betaW) W2  -0.2216027 -0.1933744
```

The dyad jackknife uses normal intervals from the jackknife standard
error, because the delete-one fits are not a bootstrap sampling
distribution:

``` r

confint(fit, boot = br_dyad)
#>                  2.5 %     97.5 %
#> (Z) Z1      -0.1995158 -0.1344302
#> (Z) Z2       0.4273622  0.5197281
#> (alphaW) W2  0.1715643  1.2976798
#> (betaW) W1  -0.6713188 -0.0287776
#> (betaW) W2  -0.3759329 -0.0295200
```

For reporting, synthesize the uncertainty checks in one table instead of
presenting four disconnected printouts. The rendered version below is a
mechanics check: the block bootstrap used only `R = 20` to keep the
vignette fast. Re-run the same table with `R = 500`–`1000` before
treating the block intervals as evidence.

``` r

ci_classical <- confint(fit, se.type = "classical")
ci_cluster <- confint(fit)
ci_block <- confint(fit, boot = br)
ci_dyad <- confint(fit, boot = br_dyad)

sign_stable <- function(ci) sign(ci[, 1]) == sign(ci[, 2])
interval_pattern <- ifelse(
    sign_stable(ci_cluster) & sign_stable(ci_dyad),
    "same sign in cluster and dyad intervals",
    ifelse(sign_stable(ci_cluster), "same sign in cluster interval only",
           "interval crosses zero")
)

data.frame(
    term = names(coef(fit)),
    estimate = round(unname(coef(fit)), 3),
    classical = paste0("[", round(ci_classical[, 1], 3), ", ", round(ci_classical[, 2], 3), "]"),
    cluster = paste0("[", round(ci_cluster[, 1], 3), ", ", round(ci_cluster[, 2], 3), "]"),
    block = paste0("[", round(ci_block[, 1], 3), ", ", round(ci_block[, 2], 3), "]"),
    dyad = paste0("[", round(ci_dyad[, 1], 3), ", ", round(ci_dyad[, 2], 3), "]"),
    cluster_over_classical = round(unname(se_cluster / se_classical), 2),
    interval_pattern = interval_pattern,
    block_note = "R = 20; mechanics only",
    row.names = NULL
)
#>          term estimate        classical          cluster            block
#> 1      (Z) Z1   -0.167 [-0.181, -0.153] [-0.179, -0.155] [-0.182, -0.156]
#> 2      (Z) Z2    0.474   [0.461, 0.486]    [0.457, 0.49]   [0.455, 0.487]
#> 3 (alphaW) W2    0.735   [0.683, 0.787]   [0.689, 0.781]   [0.677, 0.775]
#> 4  (betaW) W1   -0.350   [-0.36, -0.34] [-0.359, -0.341] [-0.361, -0.343]
#> 5  (betaW) W2   -0.203 [-0.218, -0.188] [-0.214, -0.192] [-0.222, -0.193]
#>               dyad cluster_over_classical
#> 1   [-0.2, -0.134]                   0.88
#> 2    [0.427, 0.52]                   1.26
#> 3   [0.172, 1.298]                   0.89
#> 4 [-0.671, -0.029]                   0.93
#> 5  [-0.376, -0.03]                   0.72
#>                          interval_pattern             block_note
#> 1 same sign in cluster and dyad intervals R = 20; mechanics only
#> 2 same sign in cluster and dyad intervals R = 20; mechanics only
#> 3 same sign in cluster and dyad intervals R = 20; mechanics only
#> 4 same sign in cluster and dyad intervals R = 20; mechanics only
#> 5 same sign in cluster and dyad intervals R = 20; mechanics only
```

Each bootstrap replicate refits the full model; non-converged replicates
are dropped and reported, so check the valid-replicate count. Use
`R = 500`–`1000` for publication-quality block or parametric bootstrap
intervals. The dyad jackknife ignores `R` and instead uses one refit per
deleted actor (or sender and receiver actor for rectangular fits).

## Important caveats on dependence

The classical and HC0-robust SE types assume **independence across
dyads**. Relational data violate this: the $`(i, j)`$ and $`(j, i)`$
directed dyads share members, and degree heterogeneity correlates every
tie involving the same actor. The HC0-robust SEs only correct for
heteroskedasticity/overdispersion (note how little they moved above).
The period/block bootstrap resamples whole time periods as independent
blocks: it preserves within-period network dependence, but not serial
dependence across adjacent periods.

The package ships two dyad-aware tools that *do* address it, both shown
earlier in this vignette:

- [`vcov()`](https://rdrr.io/r/stats/vcov.html) /
  [`confint()`](https://rdrr.io/r/stats/confint.html) — a multiway
  cluster-robust sandwich (sender + receiver + time), the default
  reporting interval for supported static directed-network fits with a
  stable Hessian.
- `boot_sir(type = "dyad")` — a delete-one-actor jackknife (also the
  only inference path for full-bilinear bipartite fits).

In practice:

- Report **cluster-robust** SEs
  ([`vcov()`](https://rdrr.io/r/stats/vcov.html)) for supported static
  fits with reliable analytic covariance, and use the
  bootstrap/jackknife as a cross-check.
- For dynamic `W`, unstable Hessians, or full-bilinear bipartite fits,
  prefer the dyad jackknife and response-scale sensitivity checks over
  cluster SEs.
- Do not use dyad-independent $`p`$-values as the primary reporting
  target for influence coefficients; report cluster or dyad-jackknife
  intervals and the response-scale scenario table.

## Fixed-receiver model and model comparison

Setting `fix_receiver = TRUE` constrains $`B = I`$ and estimates only
sender-side influence. This removes the alpha/beta scaling ambiguity by
dropping the receiver-side channel, but the remaining coefficients still
need adequate design rank and signal. It is faster and
better-conditioned, and is a useful descriptive baseline to compare
against the full bilinear model on the same response and mask.

``` r

fit_fr <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z,
              family = "poisson", fix_receiver = TRUE, calc_se = TRUE, seed = 1)
fr_ci <- confint(fit_fr, se.type = "cluster")
data.frame(
    term = names(coef(fit_fr)),
    estimate = round(unname(coef(fit_fr)), 3),
    cluster_low = round(fr_ci[, 1], 3),
    cluster_high = round(fr_ci[, 2], 3),
    row.names = NULL
)
#>          term estimate cluster_low cluster_high
#> 1      (Z) Z1   -0.202      -0.239       -0.164
#> 2      (Z) Z2    0.554       0.491        0.617
#> 3 (alphaW) W1   -0.192      -1.221        0.838
#> 4 (alphaW) W2   -0.260      -0.782        0.262
```

Compare specifications with fit criteria, convergence diagnostics, and a
common response-scale scenario. Here the scenario shifts the first
direct covariate by one standard deviation under each fitted model:

``` r

scenario_change <- function(model) {
    Zcf <- dat$Z
    Zcf[, , 1, ] <- Zcf[, , 1, ] + stats::sd(Zcf[, , 1, ], na.rm = TRUE)
    mu0 <- predict(model)
    mu1 <- predict(model, newdata = list(W = dat$W, X = dat$X, Z = Zcf))
    mean(mu1 - mu0, na.rm = TRUE)
}

data.frame(
    model = c("full bilinear", "fix_receiver"),
    logLik = round(c(as.numeric(logLik(fit)), as.numeric(logLik(fit_fr))), 2),
    AIC   = round(c(AIC(fit), AIC(fit_fr)), 2),
    BIC   = round(c(BIC(fit), BIC(fit_fr)), 2),
    converged = c(fit$convergence, fit_fr$convergence),
    se_reliable = c(fit$se_reliable, fit_fr$se_reliable),
    mean_scenario_change = round(c(scenario_change(fit), scenario_change(fit_fr)), 3)
)
#>           model    logLik      AIC      BIC converged se_reliable
#> 1 full bilinear -19482.18 38974.36 39012.29      TRUE        TRUE
#> 2  fix_receiver -22712.67 45433.34 45463.68      TRUE        TRUE
#>   mean_scenario_change
#> 1               -0.209
#> 2               -0.219
```

Lower AIC/BIC is better, but model choice should also ask whether the
simpler model changes the substantive response-scale conclusion.

## Out-of-sample evaluation

In-sample criteria (AIC/BIC) reward complexity; out-of-sample accuracy
is the more demanding test.
[`forecast()`](https://generics.r-lib.org/reference/forecast.html)
produces iterated plug-in forecasts, and
[`cv_sir()`](https://netify-dev.github.io/sir/reference/cv_sir.md) runs
rolling-origin (expanding-window) cross-validation, scoring each
one-step forecast against the held-out actuals alongside a last-value
naive baseline.
[`cv_sir()`](https://netify-dev.github.io/sir/reference/cv_sir.md)
conditions on the future `Z` and dynamic-`W` slices stored on the fitted
object, so treat its scores as a conditional validation target: the
future covariates must be known, pre-specified, or generated from a
scenario available at the forecast origin.

``` r

# conditional one-step forecast: fit through T-1, then pass a Z_T slice that is
# known, fixed, or scenario-generated at the forecast origin. In this simulated
# example, the held-out Z_T is used to show the API.
T_hold <- dim(dat$Y)[3] - 1
fit_train <- sir(dat$Y[, , 1:T_hold, drop = FALSE],
                 W = dat$W,
                 X = dat$X[, , 1:T_hold, drop = FALSE],
                 Z = dat$Z[, , , 1:T_hold, drop = FALSE],
                 family = "poisson", fix_receiver = TRUE,
                 calc_se = FALSE, seed = 1)
Z_hold <- dat$Z[, , , T_hold + 1, drop = FALSE]
fc <- forecast(fit_train, h = 1, Z_future = Z_hold)
model_score <- score_sir(dat$Y[, , T_hold + 1, drop = FALSE], fc, "poisson")
naive_score <- score_sir(
    dat$Y[, , T_hold + 1, drop = FALSE],
    array(dat$Y[, , T_hold], dim = dim(fc)),
    "poisson"
)
data.frame(
    metric = names(model_score),
    model = round(unname(model_score), 3),
    naive = round(unname(naive_score), 3),
    improvement = round(unname(naive_score - model_score), 3),
    row.names = NULL
)
#>     metric model  naive improvement
#> 1     rmse 1.399  2.008       0.609
#> 2      mae 0.976  1.385       0.409
#> 3 deviance 1.379 32.591      31.213

# rolling-origin conditional cross-validation: refit on 1..o, score period o + 1
cv <- cv_sir(fit_fr, initial = 60, origins = c(60, 70, 79))
data.frame(
    metric = c("rmse", "mae", "deviance"),
    model = round(unname(cv$aggregate[c("rmse", "mae", "deviance")]), 3),
    naive = round(unname(cv$aggregate[paste0("naive_", c("rmse", "mae", "deviance"))]), 3),
    improvement = round(unname(cv$aggregate[paste0("naive_", c("rmse", "mae", "deviance"))] -
                               cv$aggregate[c("rmse", "mae", "deviance")]), 3),
    effective_origins = unname(cv$eff_n[c("rmse", "mae", "deviance")]),
    row.names = NULL
)
#>     metric model  naive improvement effective_origins
#> 1     rmse 1.565  2.073       0.508                 3
#> 2      mae 1.023  1.374       0.351                 3
#> 3 deviance 1.511 31.466      29.955                 3
```

Lower RMSE/MAE/deviance scores are better. This example uses only three
origins, so it demonstrates the workflow rather than providing a final
validation study; with real applications, use more origins or a
substantively chosen rolling schedule.

[`score_sir()`](https://netify-dev.github.io/sir/reference/score_sir.md)
exposes the same family-appropriate scores (RMSE / MAE plus the Poisson
deviance) for any observed/predicted pair, so you can assemble custom
backtests.
