# Inference and model comparison

Inference for bilinear models takes more care than for an ordinary GLM.
Network data violate the independent-observations assumption behind
naive Wald standard errors — ties that share an actor are correlated —
so this vignette is mostly about getting *honest intervals* in the face
of that dependence. (A separate risk, weak identification of the
influence term, shows up as an ill-conditioned Hessian and is flagged by
the `se_reliable` field, covered below.) It walks through classical and
robust variance estimates, confidence intervals, bootstrap checks, and
information-criterion model comparison. For basic fitting see
[`vignette("sir_overview")`](https://netify-dev.github.io/sir/articles/sir_overview.md).

### Getting standard errors, in one paragraph

Standard errors are computed by default (`calc_se = TRUE`), so you do
not have to ask for them. **`summary(fit)` prints them, and
`confint(fit)`, `vcov(fit)`, and `tidy(fit)` report cluster-robust
intervals by default.** Cluster-robust is the right default for
relational data: because ties that share an actor are correlated,
classical (Hessian-based) errors that assume independent dyads can badly
overstate significance. The coefficient table in
[`summary()`](https://rdrr.io/r/base/summary.html)/[`print()`](https://rdrr.io/r/base/print.html)
shows the *classical* SE for a quick look (and labels it as such), but
the interval you should report is the cluster-robust one from
[`confint()`](https://rdrr.io/r/stats/confint.html). If the analytic
standard errors cannot be formed — a singular or ill-conditioned
Hessian, or a full-bilinear bipartite fit, which has no closed-form
covariance —
[`sir()`](https://netify-dev.github.io/sir/reference/sir.md)
**automatically falls back to the delete-one-actor jackknife**, and all
four accessors use it transparently; the fit then carries
`fit$se_source == "jackknife"` and says so in its printout. Everything
below is the *why* behind these defaults, and how to cross-check them.

| Tool | What it is useful for | Main limitation |
|:---|:---|:---|
| Classical SE | Fast model-based Wald intervals; what [`summary()`](https://rdrr.io/r/base/summary.html)/[`print()`](https://rdrr.io/r/base/print.html) show | Assumes independent dyad-period scores and a stable Hessian; can overstate significance |
| Cluster-robust SE (default) | Shared-actor dyadic dependence (each cell scored onto both actors); what [`confint()`](https://rdrr.io/r/stats/confint.html)/[`vcov()`](https://rdrr.io/r/stats/vcov.html)/[`tidy()`](https://generics.r-lib.org/reference/tidy.html) report | Still uses the Hessian as bread, so weak identification remains a problem |
| HC0 robust SE | Heteroskedasticity or overdispersion | Does not address shared-actor network dependence |
| Parametric bootstrap | Resimulates from the fitted model; confirms the model-based SEs | Assumes the model is correct, so it does not capture shared-actor dependence |
| Dyad jackknife | The automatic fallback when no analytic covariance exists (full-bilinear bipartite, ill-conditioned Hessian); also a conservative robustness check | Wider than the cluster interval; with few actors, deleting one is a large perturbation |

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

The estimator recovers these known parameters
([`vignette("sir_overview")`](https://netify-dev.github.io/sir/articles/sir_overview.md)
shows the recovery check); here we take the fit as given and focus on
putting honest *intervals* around it. In the tables below, the
coefficient labels tag each parameter’s role — `(Z)` a direct covariate
effect, `(alphaW)` a sender-influence weight, `(betaW)` a
receiver-influence weight — and the first sender weight is fixed at 1
for identifiability, so there is no `(alphaW) W1` row.

## Classical (Hessian-based) standard errors

`sir(calc_se = TRUE)` computes classical standard errors from the
inverse observed-information matrix.
[`summary()`](https://rdrr.io/r/base/summary.html) prints them for a
quick look, along with a note that they are Hessian-based and can
overstate significance under dyadic dependence:

``` r

summary(fit)
#> 
#> ── Social Influence Regression Model ───────────────────────────────────────────
#> Network: 14 nodes, 80 time periods (directed)
#> Family: "poisson" | Method: "ALS"
#> Observations: 14560
#> ────────────────────────────────────────────────────────────────────────────────
#> 
#> ── Coefficients ──
#> 
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Z) Z1      -0.166973   0.007030  -23.75   <2e-16 ***
#> (Z) Z2       0.473545   0.006552   72.27   <2e-16 ***
#> (alphaW) W2  0.734622   0.026506   27.71   <2e-16 ***
#> (betaW) W1  -0.350048   0.005033  -69.55   <2e-16 ***
#> (betaW) W2  -0.202726   0.007608  -26.65   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Std. Error, z value, Pr(>|z|) and significance stars above are classical
#> (Hessian-based) and may overstate significance under dyadic dependence. NOTE:
#> `confint(fit)`/`vcov(fit)`/`tidy(fit)` default to cluster-robust inference
#> (actor-clustered, t(G-1) reference) --- the same estimator for directed and
#> symmetric fits; request classical intervals with `confint(se.type =
#> "classical")` if dyads are independent.
#> (Z) = direct covariate effect; (alphaW) = sender influence; (betaW) = receiver
#> influence
#> ────────────────────────────────────────────────────────────────────────────────
#> 
#> ── Model Fit ──
#> 
#> • Log-Likelihood: -19482.18
#> • AIC: 38974.36
#> • BIC: 39012.29
#> • Residual deviance: 15933.43
#> • Dispersion (Pearson chi-sq / resid df): 1.022
#> ✔ Converged in 3 iterations
#> ────────────────────────────────────────────────────────────────────────────────
#> 
#> ── Influence Matrices ──
#> 
#> A matrix (sender effects):
#> • Mean: 0.0696
#> • SD: 1.2133
#> • Range: [-3.6025, 3.0012]
#> B matrix (receiver effects):
#> • Mean: -0.0278
#> • SD: 0.3912
#> • Range: [-1.0295, 1.1373]
```

The closing **Influence Matrices** block summarizes the fitted
$`m \times m`$ sender and receiver matrices $`A`$ and $`B`$ (each a
weighted sum of the `W` slices), so their entries live on a wider scale
than the coefficient table above;
[`vignette("sir_overview")`](https://netify-dev.github.io/sir/articles/sir_overview.md)
covers how to read them.

For reporting, the default accessors use cluster-robust uncertainty when
the fit supports it; request the classical SEs explicitly when you need
a Hessian-only comparison:

``` r

round(sqrt(diag(vcov(fit, type = "classical"))), 4)   # the SEs reported by summary()
#>      (Z) Z1      (Z) Z2 (alphaW) W2  (betaW) W1  (betaW) W2 
#>      0.0070      0.0066      0.0265      0.0050      0.0076
confint(fit, se.type = "classical")                   # Hessian-only Wald intervals
#>                  2.5 %     97.5 %
#> (Z) Z1      -0.1807515 -0.1531944
#> (Z) Z2       0.4607029  0.4863873
#> (alphaW) W2  0.6826715  0.7865726
#> (betaW) W1  -0.3599133 -0.3401830
#> (betaW) W2  -0.2176376 -0.1878153
```

## Robust (sandwich) standard errors

The HC0 sandwich estimator $`H^{-1} S H^{-1}`$ (with $`S`$ the empirical
score covariance) corrects for heteroskedasticity or overdispersion. Our
simulation is a correctly specified Poisson, so there is nothing to
correct and the robust SEs match the classical ones almost exactly — a
useful confirmation that the model is well specified:

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
```

To see HC0 do real work, refit on overdispersed counts — here redrawn
from a negative binomial, whose variance far exceeds its mean. Now the
robust SEs come out about three times the classical ones, correctly
widening for the extra dispersion the Poisson likelihood would otherwise
ignore:

``` r

set.seed(99)
Y_od <- dat$Y
Y_od[] <- rnbinom(length(Y_od), mu = pmax(dat$Y, 0.1), size = 0.5)   # heavier-than-Poisson counts
fit_od <- sir(Y_od, W = dat$W, X = dat$X, Z = dat$Z,
              family = "poisson", calc_se = TRUE, seed = 1)

data.frame(
    term = names(coef(fit_od)),
    robust_over_classical = round(sqrt(diag(vcov(fit_od, type = "robust"))) /
                                  sqrt(diag(vcov(fit_od, type = "classical"))), 2),
    row.names = NULL
)
#>          term robust_over_classical
#> 1      (Z) Z1                  3.05
#> 2      (Z) Z2                  3.18
#> 3 (alphaW) W2                  3.06
#> 4  (betaW) W1                  3.19
#> 5  (betaW) W2                  3.11
```

HC0 handles dispersion and heteroskedasticity, but neither HC0 form
models the shared-actor dependence that the **cluster-robust** estimator
addresses next.

## Cluster-robust standard errors

Every tie involving an actor is correlated, so
[`vcov()`](https://rdrr.io/r/stats/vcov.html) clusters on the
**actors**: each dyad-period score is assigned to both of its endpoint
actors and summed within actor. (Mechanically, an HC1 small-sample
factor and a $`t(G - 1)`$ reference, with $`G`$ the number of actors,
complete the sandwich.) This is the default reporting interval for
supported static directed and symmetric fits when the Hessian bread is
stable; report it in place of classical Hessian SEs:

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
#> 1      (Z) Z1    0.0070  0.0101  1.44
#> 2      (Z) Z2    0.0066  0.0123  1.88
#> 3 (alphaW) W2    0.0265  0.0277  1.04
#> 4  (betaW) W1    0.0050  0.0072  1.42
#> 5  (betaW) W2    0.0076  0.0091  1.19
```

Here the cluster-robust SEs are 1.0–1.9 times the classical ones:
because they stop assuming independent dyad-times, they are the more
conservative — and the ones to report. `confint(fit)` returns the
matching Wald intervals:

``` r

confint(fit)
#>                  2.5 %     97.5 %
#> (Z) Z1      -0.1887952 -0.1451508
#> (Z) Z2       0.4469773  0.5001130
#> (alphaW) W2  0.6748119  0.7944323
#> (betaW) W1  -0.3655139 -0.3345825
#> (betaW) W2  -0.2223161 -0.1831368
```

On this simulation the gap is modest because the periods are nearly
independent. On real relational data it is not. Refitting on the bundled
`icews` inter-state conflict network, the cluster SEs come out many
times the classical ones:

``` r

data(icews)
ic <- 1:10; ip <- 1:24
ifit <- sir(icews$Y[ic, ic, ip, drop = FALSE], W = icews$W[ic, ic, , drop = FALSE],
            X = icews$X[ic, ic, ip, drop = FALSE], Z = icews$Z[ic, ic, , ip, drop = FALSE],
            family = "poisson", seed = 1, max_iter = 20)

data.frame(
    term = names(coef(ifit)),
    cluster_over_classical = round(sqrt(diag(vcov(ifit, type = "cluster"))) /
                                   sqrt(diag(vcov(ifit, type = "classical"))), 1),
    row.names = NULL
)
#>                   term cluster_over_classical
#> 1            (Z) mConf                   18.3
#> 2         (Z) mConf_ji                   21.8
#> 3       (Z) minDistLog                   43.3
#> 4             (Z) ally                   23.2
#> 5         (Z) verbCoop                   80.7
#> 6        (alphaW) ally                   34.1
#> 7    (alphaW) verbCoop                   21.0
#> 8  (alphaW) minDistLog                   27.0
#> 9          (betaW) int                   49.1
#> 10        (betaW) ally                   24.5
#> 11    (betaW) verbCoop                   66.5
#> 12  (betaW) minDistLog                   12.4
```

Real conflict data carry exactly the shared-actor dependence the
classical SEs assume away, so the classical intervals would badly
overstate significance — the cluster SEs are an order of magnitude
larger. This is why the cluster interval is the default and the reported
one.

## Bootstrap inference

[`boot_sir()`](https://netify-dev.github.io/sir/reference/boot_sir.md)
refits the model on resampled data and forms distribution-free intervals
that do not assume the Wald normal shape. Two are worth running: the
**parametric** bootstrap resimulates outcomes from the fitted model (a
model-based cross-check), and the **dyad jackknife** deletes each actor
in turn (a conservative robustness check, and the fallback when analytic
SEs are unavailable).

``` r

br_par  <- boot_sir(fit, R = 200, type = "parametric", seed = 123, trace = FALSE)
br_dyad <- boot_sir(fit, type = "dyad", seed = 123, trace = FALSE)

data.frame(
    method = c("Parametric bootstrap", "Dyad jackknife"),
    valid_refits = c(br_par$n_valid, br_dyad$n_valid),
    total_refits = c(br_par$n_total, br_dyad$n_total),
    interval_type = c(br_par$interval, br_dyad$interval)
)
#>                 method valid_refits total_refits    interval_type
#> 1 Parametric bootstrap          200          200       percentile
#> 2       Dyad jackknife           14           14 normal-jackknife
```

Each replicate refits the full model and non-converged ones are dropped,
so check the valid-replicate count. The **parametric** bootstrap
intervals track the classical ones and every one excludes zero — a
distribution-free confirmation that the estimates are well separated
from zero:

``` r

confint(fit, boot = br_par)
#>                  2.5 %     97.5 %
#> (Z) Z1      -0.1809816 -0.1565513
#> (Z) Z2       0.4624189  0.4887205
#> (alphaW) W2  0.6774419  0.7810377
#> (betaW) W1  -0.3590865 -0.3390764
#> (betaW) W2  -0.2167343 -0.1874963
```

The **dyad jackknife** is deliberately conservative: deleting a whole
actor is a large perturbation (especially with only `m = 14` actors), so
its intervals run several times wider, most of all for the bilinear
influence weights. Read it as a worst-case bound, not a competing
precise estimate:

``` r

confint(fit, boot = br_dyad)
#>                  2.5 %     97.5 %
#> (Z) Z1      -0.1995158 -0.1344302
#> (Z) Z2       0.4273622  0.5197281
#> (alphaW) W2  0.1715643  1.2976798
#> (betaW) W1  -0.6713188 -0.0287776
#> (betaW) W2  -0.3759329 -0.0295200
```

For reporting, place the intervals side by side — the classical (which
ignores dyadic dependence), the cluster-robust default, and the
conservative dyad jackknife:

``` r

fmt <- function(ci) sprintf("[%.3f, %.3f]", ci[, 1], ci[, 2])

data.frame(
    term           = names(coef(fit)),
    estimate       = round(unname(coef(fit)), 3),
    classical      = fmt(confint(fit, se.type = "classical")),
    cluster        = fmt(confint(fit)),
    dyad_jackknife = fmt(confint(fit, boot = br_dyad)),
    row.names = NULL
)
#>          term estimate        classical          cluster   dyad_jackknife
#> 1      (Z) Z1   -0.167 [-0.181, -0.153] [-0.189, -0.145] [-0.200, -0.134]
#> 2      (Z) Z2    0.474   [0.461, 0.486]   [0.447, 0.500]   [0.427, 0.520]
#> 3 (alphaW) W2    0.735   [0.683, 0.787]   [0.675, 0.794]   [0.172, 1.298]
#> 4  (betaW) W1   -0.350 [-0.360, -0.340] [-0.366, -0.335] [-0.671, -0.029]
#> 5  (betaW) W2   -0.203 [-0.218, -0.188] [-0.222, -0.183] [-0.376, -0.030]
```

The estimate is shared across the columns (one fit); what changes is the
interval. Classical is the narrowest (it assumes independent dyads), the
cluster default widens it for shared-actor dependence, and the jackknife
is the conservative outer bound. Report the cluster interval and use the
jackknife as a robustness check. Use `R = 500`–`1000` for
publication-quality bootstrap intervals.

## Choosing an approach

The choice reduces to three rules:

- **Report cluster-robust intervals**
  ([`vcov()`](https://rdrr.io/r/stats/vcov.html) /
  [`confint()`](https://rdrr.io/r/stats/confint.html), the default) for
  directed, symmetric, and dynamic (4D) `W` fits with a stable Hessian
  (`se_reliable = TRUE`).
- **Let the jackknife fallback handle the hard cases.** When no analytic
  covariance exists — full-bilinear bipartite fits, or an
  ill-conditioned Hessian —
  [`sir()`](https://netify-dev.github.io/sir/reference/sir.md) already
  attaches the delete-one-actor jackknife automatically, so
  [`confint()`](https://rdrr.io/r/stats/confint.html)/[`vcov()`](https://rdrr.io/r/stats/vcov.html)/[`tidy()`](https://generics.r-lib.org/reference/tidy.html)
  keep working and `fit$se_source` reads `"jackknife"`; you do not need
  to call
  [`boot_sir()`](https://netify-dev.github.io/sir/reference/boot_sir.md)
  yourself. (You still can, e.g. `boot_sir(type = "dyad")` then
  `confint(fit, boot = bs)`, if you want to inspect the resampling.)
- **Do not** report dyad-independent classical $`p`$-values as the
  headline result for influence coefficients; pair the cluster interval
  with a response-scale scenario (below).

### How to tell which standard errors you have

`fit$se_source` is `"jackknife"` when the automatic fallback fired, and
absent (NULL) otherwise — meaning analytic SEs are available and
[`confint()`](https://rdrr.io/r/stats/confint.html)/[`vcov()`](https://rdrr.io/r/stats/vcov.html)/[`tidy()`](https://generics.r-lib.org/reference/tidy.html)
report the cluster-robust sandwich (including for dynamic, 4D `W`). The
`summary(fit)` and `print(fit)` printouts also label the standard errors
they show (classical in the table, jackknife when the fallback fired),
so you are never left guessing which kind you are looking at.

## Fixed-receiver model and model comparison

Setting `fix_receiver = TRUE` constrains $`B = I`$ and estimates only
sender-side influence. It is faster and better-conditioned, and serves
as a constrained baseline. Because our simulation has a genuine
receiver-side channel that this constraint drops, we expect the full
model to fit better — so the useful question is not what the constrained
model’s influence coefficients are, but whether dropping the receiver
side costs fit. We answer that with information criteria rather than by
reading the constrained coefficients:

``` r

fit_fr <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z,
              family = "poisson", fix_receiver = TRUE, calc_se = TRUE, seed = 1)
```

Compare specifications with fit criteria and a common response-scale
scenario. Because the two models differ in how the *lagged network*
propagates, we shift the lagged signal `X` up by one standard deviation
and read how each model’s predicted mean responds:

``` r

scenario_change <- function(model) {
    Xcf <- dat$X + stats::sd(dat$X, na.rm = TRUE)
    mu0 <- predict(model)
    mu1 <- predict(model, newdata = list(W = dat$W, X = Xcf, Z = dat$Z))
    mean(mu1 - mu0, na.rm = TRUE)
}

data.frame(
    model = c("full bilinear", "fix_receiver"),
    logLik = round(c(as.numeric(logLik(fit)), as.numeric(logLik(fit_fr))), 2),
    AIC   = round(c(AIC(fit), AIC(fit_fr)), 2),
    BIC   = round(c(BIC(fit), BIC(fit_fr)), 2),
    se_reliable = c(fit$se_reliable, fit_fr$se_reliable),
    mean_response_to_X = round(c(scenario_change(fit), scenario_change(fit_fr)), 3)
)
#>           model    logLik      AIC      BIC se_reliable mean_response_to_X
#> 1 full bilinear -19482.18 38974.36 39012.29        TRUE              0.650
#> 2  fix_receiver -22712.67 45433.34 45463.68        TRUE              0.001
```

The `se_reliable` column flags whether the Hessian was well-conditioned
enough to trust the analytic SEs (both fits here are reliable). Lower
AIC/BIC is better, and the full bilinear model wins decisively —
confirming the receiver-side channel carries real signal. The scenario
makes the difference concrete: shifting the lagged network moves the
full model’s predicted mean noticeably, but barely moves the
fixed-receiver model’s. Forcing $`B = I`$ leaves a sender map (the `A`
matrix) that is both weaker in spread and nearly self-cancelling, so a
uniform shift in the lagged network produces almost no net change in its
predicted mean, whereas the full model’s stronger, coherent influence
structure responds. The much larger drop in log-likelihood for the
constrained model confirms the receiver side carries signal it cannot
recover.
[`vignette("sir_overview")`](https://netify-dev.github.io/sir/articles/sir_overview.md)
covers `predict(newdata = )` scenarios in full.

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
                 family = "poisson", calc_se = FALSE, seed = 1)
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
#> 1     rmse 1.130  2.008       0.878
#> 2      mae 0.846  1.385       0.539
#> 3 deviance 1.077 32.591      31.515

# rolling-origin conditional cross-validation: initial is the minimum training
# length and origins are the cut-points o (each refits on 1..o and scores o + 1)
cv <- cv_sir(fit, initial = 60, origins = c(60, 70, 79))
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
#> 1     rmse 1.146  2.073       0.927                 3
#> 2      mae 0.867  1.374       0.506                 3
#> 3 deviance 1.139 31.466      30.327                 3
```

Lower RMSE/MAE/deviance scores are better. The model improves all three
over the naive last-value baseline; the deviance gap is the largest
because Poisson deviance penalizes the naive predictor heavily when its
carried-forward counts miss, while RMSE and MAE are gentler. This
example uses only three origins, so it demonstrates the workflow rather
than providing a final validation study; with real applications, use
more origins or a substantively chosen rolling schedule.

[`score_sir()`](https://netify-dev.github.io/sir/reference/score_sir.md)
exposes the same family-appropriate scores (RMSE / MAE plus the Poisson
deviance) for any observed/predicted pair, so you can assemble custom
backtests.
