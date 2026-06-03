# Distribution families, network types, and data preparation

The SIR framework handles a range of outcome types, network structures,
and covariate configurations beyond the Poisson directed-network case in
the overview. This vignette covers the distribution families, symmetric
and bipartite networks, dynamic (time-varying) influence covariates, and
the data-preparation utilities. For basic fitting see
[`vignette("sir_overview")`](https://netify-dev.github.io/sir/articles/sir_overview.md).

## Distribution families

The distribution-family examples below use
[`sim_sir()`](https://netify-dev.github.io/sir/reference/sim_sir.md), so
each fits data with a *known* influence signal and the coefficients are
interpretable. The later sections (symmetric, bipartite, dynamic) use
small illustrative networks to show the *call mechanics* of each
structure; there the simulated outcome is plain noise, so the fitted
coefficients are not meaningful.

### Normal (continuous outcomes)

For continuous relational data (trade volumes, sentiment scores) the
Normal family uses an identity link. The influence parameters are
interpreted exactly as in the Poisson case: $`\alpha`$ and $`\beta`$
identify which covariates account for sender and receiver influence,
with effects entering the linear predictor.

``` r

dn  <- sim_sir(m = 14, T_len = 80, p = 2, q = 2, family = "normal", seed = 11)
fit_norm <- sir(dn$Y, W = dn$W, X = dn$X, Z = dn$Z,
                family = "normal", calc_se = FALSE, seed = 1)
target_norm <- c(dn$theta, dn$alpha[-1], dn$beta)
data.frame(
    term = names(coef(fit_norm)),
    estimate = round(unname(coef(fit_norm)), 3),
    target = round(target_norm, 3),
    abs_error = round(abs(unname(coef(fit_norm)) - target_norm), 3),
    row.names = NULL
)
#>          term estimate target abs_error
#> 1      (Z) Z1   -0.411 -0.414     0.002
#> 2      (Z) Z2   -0.207 -0.210     0.003
#> 3 (alphaW) W2   -0.156 -0.177     0.022
#> 4  (betaW) W1    0.000  0.008     0.008
#> 5  (betaW) W2   -0.461 -0.455     0.006
```

The synthetic data were generated from a known Normal SIR process, so
the coefficient signs and magnitudes are meaningful for this example. We
show recovery against the generating values rather than a star table
because this section is about family behavior, not inference.

On the response scale, the Normal family returns fitted means directly:

``` r

mu_norm <- predict(fit_norm)
round(quantile(mu_norm[!is.na(mu_norm)], c(0.05, 0.5, 0.95)), 3)
#>     5%    50%    95% 
#> -1.177 -0.028  1.138
```

### Binomial (binary outcomes)

For binary relational outcomes (presence/absence of a tie, conflict
onset) the Binomial family uses a logit link. The full bilinear
contribution $`(A X_t B^\top)_{ij}`$ enters the log-odds, so an
individual $`\alpha`$ or $`\beta`$ coefficient is not a standalone
odds-ratio term; its effect depends on `X`, the other influence side,
and the covariate scale.

``` r

db  <- sim_sir(m = 16, T_len = 120, p = 2, q = 1, family = "binomial", seed = 12)
fit_bin <- sir(db$Y, W = db$W, X = db$X, Z = db$Z,
               family = "binomial", calc_se = FALSE, seed = 1)
target_bin <- c(db$theta, db$alpha[-1], db$beta)
data.frame(
    term = names(coef(fit_bin)),
    estimate = round(unname(coef(fit_bin)), 3),
    target = round(target_bin, 3),
    abs_error = round(abs(unname(coef(fit_bin)) - target_bin), 3),
    row.names = NULL
)
#>          term estimate target abs_error
#> 1      (Z) Z1   -0.304 -0.321     0.018
#> 2 (alphaW) W2   -0.484 -0.444     0.040
#> 3  (betaW) W1    0.488  0.473     0.014
#> 4  (betaW) W2   -0.272 -0.287     0.015
```

Again, these are simulated data with a known generating process, so
recovery and response-scale probabilities are the useful outputs.

[`predict()`](https://rdrr.io/r/stats/predict.html) returns fitted tie
probabilities by default. Use `predict(type = "link")` and then
[`plogis()`](https://rdrr.io/r/stats/Logistic.html) only when you
explicitly request the log-odds scale:

``` r

mu_bin <- predict(fit_bin)                 # fitted tie probabilities
round(quantile(mu_bin[!is.na(mu_bin)], c(0.05, 0.5, 0.95)), 3)
#>    5%   50%   95% 
#> 0.290 0.496 0.703
```

## Symmetric (undirected) networks

For undirected networks where $`y_{i,j} = y_{j,i}`$, set
`symmetric = TRUE`. The model fits the upper triangle and constrains
$`B = I`$. This is an upper-triangle, sender-side representation of
undirected data, not a fully order-invariant undirected bilinear model.
For Poisson and Binomial outcomes, provide a symmetric `Y` directly;
averaging asymmetric discrete matrices can create invalid non-integer or
non-binary observations. The influence covariates in `W` must also be
symmetric and meaningful for unordered pairs. This suits undirected
trade, alliance, or co-sponsorship networks when that representation is
substantively acceptable.

``` r

set.seed(1)
m <- 14; T_len <- 60; p <- 2
W <- array(0, dim = c(m, m, p))
for (k in seq_len(p)) {
    W_k <- matrix(rnorm(m * m), m, m)
    W[, , k] <- (W_k + t(W_k)) / 2
    diag(W[, , k]) <- 0
}

Y_sym <- array(0, dim = c(m, m, T_len))
for (t in 1:T_len) {
    Y_t <- matrix(0, m, m)
    upper <- upper.tri(Y_t)
    Y_t[upper] <- rpois(sum(upper), 2)
    Y_t <- Y_t + t(Y_t)
    diag(Y_t) <- NA
    Y_sym[, , t] <- Y_t
}
X_sym <- array(0, dim = c(m, m, T_len))
for (t in 2:T_len) {
    X_sym[, , t] <- log(Y_sym[, , t - 1] + 1)
    X_sym[, , t][is.na(X_sym[, , t])] <- 0
}

fit_sym <- sir(Y_sym, W = W, X = X_sym, family = "poisson",
               symmetric = TRUE, calc_se = FALSE, seed = 1)
data.frame(
    converged = fit_sym$convergence,
    symmetric = fit_sym$symmetric,
    fix_receiver = fit_sym$fix_receiver
)
#>   converged symmetric fix_receiver
#> 1      TRUE      TRUE         TRUE
```

This chunk demonstrates the API and structural checks. The outcome was
generated as simple noise, so the fitted coefficient values are not
substantively interpretable.

> These illustrative chunks build `X` as a plain `log(Y_prev + 1)` for
> brevity. When you intend to forecast, build `X` the way
> [`sim_sir()`](https://netify-dev.github.io/sir/reference/sim_sir.md)
> and [`forecast()`](https://generics.r-lib.org/reference/forecast.html)
> assume — divide the lag by `max(m - 1, 1)` (the number of partners
> summed over) — or pass a matching `infl_scale` to
> [`forecast()`](https://generics.r-lib.org/reference/forecast.html).

## Bipartite networks

When senders and receivers are **distinct populations** (a rectangular
$`Y`$, e.g. countries $`\to`$ NGOs, or legislators co-sponsoring others’
bills), bipartite structure is detected automatically and
`fix_receiver = TRUE` is enforced.

The detail that trips people up: **influence flows on the sender side,
so `W` describes the senders.** For an $`n_1 \times n_2`$ outcome, `W`
must be $`n_1 \times n_1 \times p`$ (sender-by-sender), *not* keyed to
the receivers. `X` matches `Y` at $`n_1 \times n_2`$.

``` r

set.seed(2)
n1 <- 10; n2 <- 15; T_len <- 40; p <- 2
Y_bp <- array(rpois(n1 * n2 * T_len, 2), dim = c(n1, n2, T_len))
W_bp <- array(rnorm(n1 * n1 * p), dim = c(n1, n1, p))   # sender-by-sender
X_bp <- array(0, dim = c(n1, n2, T_len))
for (t in 2:T_len) X_bp[, , t] <- log(Y_bp[, , t - 1] + 1)
Z_bp <- array(rnorm(n1 * n2 * 1 * T_len), dim = c(n1, n2, 1, T_len))

fit_bp <- sir(Y_bp, W = W_bp, X = X_bp, Z = Z_bp,
              family = "poisson", fix_receiver = TRUE, calc_se = FALSE, seed = 1)
data.frame(
    converged = fit_bp$convergence,
    bipartite = fit_bp$bipartite,
    fix_receiver = fit_bp$fix_receiver,
    senders = fit_bp$n1,
    receivers = fit_bp$n2
)
#>   converged bipartite fix_receiver senders receivers
#> 1      TRUE      TRUE         TRUE      10        15
```

Again, this is an API demonstration: the simulated outcome is plain
noise, so only the structural output should be read.

For a *square* array whose rows and columns are nonetheless distinct
populations, set `bipartite = TRUE` explicitly so the receiver side is
not estimated.

### Full bilinear bipartite (with `W_recv`)

The default above collapses `B` to the identity. To estimate the full
bilinear model $`A X B'`$ for two-mode data — with a *separate*
receiver-side influence structure — supply `W_recv`, an
$`n_2 \times n_2 \times p_2`$ array. The sender influence `A` is then
built from `W` and the receiver influence `B` from `W_recv`, fit by
alternating GLM with `alpha_1 = 1` pinning the scale.

Here we simulate a two-mode network *with* genuine bilinear structure so
the fit has something to recover (the earlier `Y_bp` was pure noise,
useful only to show the API). True sender weights are `alpha = (1, 0.6)`
and receiver weights `beta = (0.8, -0.5)`:

``` r

set.seed(204)
n1 <- 10; n2 <- 6; Tn <- 30
Wf  <- array(rnorm(n1 * n1 * 2), dim = c(n1, n1, 2))   # sender-by-sender
Wr  <- array(rnorm(n2 * n2 * 2), dim = c(n2, n2, 2))   # receiver-by-receiver
Xf  <- array(rnorm(n1 * n2 * Tn) / sqrt(n2), dim = c(n1, n2, Tn))
A_t <- Wf[, , 1] + 0.6 * Wf[, , 2]
B_t <- 0.8 * Wr[, , 1] - 0.5 * Wr[, , 2]
Yf  <- array(0, dim = c(n1, n2, Tn))
for (t in 1:Tn) Yf[, , t] <- A_t %*% Xf[, , t] %*% t(B_t) +
    matrix(rnorm(n1 * n2, sd = 0.3), n1, n2)

fit_full <- sir(Yf, W = Wf, X = Xf, W_recv = Wr,
                family = "normal", calc_se = FALSE, seed = 1)
data.frame(
    full_bilinear = isTRUE(fit_full$full_bilinear),
    converged = fit_full$convergence
)
#>   full_bilinear converged
#> 1          TRUE      TRUE

data.frame(
    term = names(coef(fit_full)),
    estimate = round(unname(coef(fit_full)), 3),
    target = c(0.6, 0.8, -0.5),
    row.names = NULL
)
#>           term estimate target
#> 1  (alphaW) W2    0.600    0.6
#> 2 (betaWr) Wr1    0.805    0.8
#> 3 (betaWr) Wr2   -0.502   -0.5
```

Analytic standard errors are not available on this path; use the
delete-one-actor dyad jackknife for inference:

``` r

bj <- boot_sir(fit_full, type = "dyad", seed = 1)
confint(bj)
#>                   2.5 %     97.5 %
#> (alphaW) W2   0.5010128  0.6984382
#> (betaWr) Wr1  0.6968568  0.9125953
#> (betaWr) Wr2 -0.5841368 -0.4197630
```

## Dynamic influence covariates

Often the factors that mediate influence change over time — alliances
evolve, trade shifts. Supply `W` as a 4D array
($`m \times m \times p \times T`$). The coefficients
$`\boldsymbol{\alpha}`$ and $`\boldsymbol{\beta}`$ are still estimated
as time-invariant, but the reconstructed influence matrices $`A_t`$ now
vary with $`t`$ because the covariates do.

``` r

set.seed(3)
m <- 14; T_len <- 40; p <- 2
Y_dyn <- array(rpois(m * m * T_len, 2), dim = c(m, m, T_len))
for (t in 1:T_len) diag(Y_dyn[, , t]) <- NA
X_dyn <- array(0, dim = c(m, m, T_len))
for (t in 2:T_len) X_dyn[, , t] <- log(Y_dyn[, , t - 1] + 1)
X_dyn[is.na(X_dyn)] <- 0
W_dyn <- array(rnorm(m * m * p * T_len), dim = c(m, m, p, T_len))   # 4D, time-varying

fit_dyn <- sir(Y_dyn, W = W_dyn, X = X_dyn, family = "poisson",
               fix_receiver = TRUE, calc_se = FALSE, max_iter = 10, seed = 1)

dim(fit_dyn$A)        # A is now m x m x T
#> [1] 14 14 40
fit_dyn$dynamic_W
#> [1] TRUE
```

The coefficients are fixed across time, but the realized channels change
with `W_t`. Because this example uses simulated noise to demonstrate the
API, the right rendered output is a structural check rather than a
ranked channel table:

``` r

dynamic_period_check <- function(period) {
    A_t <- fit_dyn$A[, , period]
    W_t <- W_dyn[, , , period]
    diag(A_t) <- NA
    data.frame(
        period = period,
        A_dimensions = paste(dim(A_t), collapse = " x "),
        finite_off_diagonal = sum(is.finite(A_t)),
        mean_abs_A = round(mean(abs(A_t), na.rm = TRUE), 3),
        mean_abs_W = round(mean(abs(W_t), na.rm = TRUE), 3),
        row.names = NULL
    )
}
rbind(
    dynamic_period_check(1),
    dynamic_period_check(T_len)
)
#>   period A_dimensions finite_off_diagonal mean_abs_A mean_abs_W
#> 1      1      14 x 14                 182      0.014      0.834
#> 2     40      14 x 14                 182      0.013      0.798
```

## Data preparation

### Converting edge lists to arrays

Network data often arrives as a long-format edge list (sender, receiver,
time, value).
[`cast_array()`](https://netify-dev.github.io/sir/reference/cast_array.md)
converts it to the arrays
[`sir()`](https://netify-dev.github.io/sir/reference/sir.md) expects.

``` r

set.seed(4)
edge_list <- expand.grid(i = paste0("n", 1:5), j = paste0("n", 1:5), t = 1:3)
edge_list <- edge_list[edge_list$i != edge_list$j, ]
edge_list$conflict <- rpois(nrow(edge_list), lambda = 2)

Y_from_el <- cast_array(edge_list, var = "conflict")
dim(Y_from_el)
#> [1] 5 5 3
dimnames(Y_from_el)[[1]]
#> [1] "n1" "n2" "n3" "n4" "n5"
```

### Constructing relational covariates

[`rel_covar()`](https://netify-dev.github.io/sir/reference/rel_covar.md)
builds main ($`z_{ij}`$), reciprocal ($`z_{ji}`$), and transitive
($`\sum_k s_{ik} s_{kj}`$ on the symmetrized network) effects from a
base dyadic variable. These capture common relational mechanisms and
drop straight into the exogenous covariate array $`Z`$.

``` r

set.seed(5)
trade   <- array(abs(rnorm(8 * 8 * 4)), dim = c(8, 8, 4))   # base dyadic variable
Z_trade <- rel_covar(trade, "trade")
dim(Z_trade)
#> [1] 8 8 3 4
dimnames(Z_trade)[[3]]
#> [1] "trade"       "trade_recip" "trade_trans"
```

### Simulating SIR data

[`sim_sir()`](https://netify-dev.github.io/sir/reference/sim_sir.md)
(used throughout these vignettes) generates synthetic data from the SIR
process — useful for simulation studies, power analysis, and verifying
parameter recovery under controlled conditions.

``` r

dat <- sim_sir(m = 10, T_len = 8, p = 2, q = 1, family = "poisson", seed = 42)
str(dat[c("Y", "W", "X", "Z", "alpha", "beta", "theta")], max.level = 1)
#> List of 7
#>  $ Y    : num [1:10, 1:10, 1:8] 0 1 1 1 0 2 2 0 3 1 ...
#>  $ W    : num [1:10, 1:10, 1:2] 0 0.539 0.58 -0.658 1.555 ...
#>  $ X    : num [1:10, 1:10, 1:8] 0 0 0 0 0 0 0 0 0 0 ...
#>  $ Z    : num [1:10, 1:10, 1, 1:8] -0.628 -0.29 0.206 0.588 -1.024 ...
#>  $ alpha: num [1:2] 1 0.411
#>  $ beta : num [1:2] -0.169 0.109
#>  $ theta: num 0.237
```
