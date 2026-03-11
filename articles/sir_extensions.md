# Distribution families, network types, and data preparation

The SIR framework accommodates a range of outcome types, network
structures, and covariate configurations beyond the Poisson
directed-network case presented in the overview. This vignette covers
the available distribution families, symmetric and bipartite network
structures, dynamic (time-varying) influence covariates, and data
preparation utilities for converting common data formats into the arrays
required by
[`sir()`](https://netify-dev.github.io/sir/reference/sir.md). For basic
model fitting see
[`vignette("sir_overview")`](https://netify-dev.github.io/sir/articles/sir_overview.md).

## Setup

We reuse the influence covariates and exogenous covariate from the
overview.

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
```

## Distribution families

### Normal (continuous outcomes)

For continuous relational data such as trade volumes or sentiment
scores, the Normal family uses an identity link. The interpretation of
the influence parameters is the same as in the Poisson case: the
$\alpha$ and $\beta$ coefficients identify which covariates account for
sender and receiver influence patterns.

``` r
Y_norm = array(rnorm(m * m * T_len, mean = 3, sd = 1), dim = c(m, m, T_len))
for (t in 1:T_len) diag(Y_norm[,,t]) = NA

X_norm = array(0, dim = c(m, m, T_len))
for (t in 2:T_len) X_norm[,,t] = Y_norm[,,t-1]
X_norm[is.na(X_norm)] = 0

fit_norm = sir(
    Y = Y_norm, W = W, X = X_norm, Z = Z,
    family = "normal", method = "ALS", calc_se = TRUE
)
summary(fit_norm)
#>                       Estimate Std. Error z value Pr(>|z|)    
#> (Z) distance          0.016757   0.043840   0.382  0.70228    
#> (alphaW) shared_group 0.407559   0.138459   2.944  0.00324 ** 
#> (betaW) proximity     0.013482   0.001354   9.957  < 2e-16 ***
#> (betaW) shared_group  0.005168   0.001602   3.225  0.00126 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Binomial (binary outcomes)

For binary relational outcomes (e.g., presence or absence of a
diplomatic tie, conflict onset), the Binomial family uses a logit link.
The influence parameters now operate on the log-odds scale, so a
positive $\alpha$ on a covariate indicates that the corresponding
relationship characteristic increases the log-odds of influence.

``` r
Y_bin = array(rbinom(m * m * T_len, 1, 0.3), dim = c(m, m, T_len))
for (t in 1:T_len) diag(Y_bin[,,t]) = NA

X_bin = array(0, dim = c(m, m, T_len))
for (t in 2:T_len) X_bin[,,t] = Y_bin[,,t-1]
X_bin[is.na(X_bin)] = 0

fit_bin = sir(
    Y = Y_bin, W = W, X = X_bin, Z = Z,
    family = "binomial", method = "ALS", calc_se = TRUE
)
summary(fit_bin)
#>                        Estimate Std. Error z value Pr(>|z|)  
#> (Z) distance           0.008763   0.066122   0.133    0.895  
#> (alphaW) shared_group  1.769791   1.459395   1.213    0.225  
#> (betaW) proximity     -0.034493   0.016270  -2.120    0.034 *
#> (betaW) shared_group   0.006164   0.010838   0.569    0.570  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## Symmetric (undirected) networks

For undirected networks where $y_{i,j} = y_{j,i}$, setting
`symmetric = TRUE` symmetrizes $\mathbf{Y}$, fits only the upper
triangle, and constrains $\mathbf{B} = \mathbf{I}$. The sender/receiver
distinction is meaningless for undirected data, so the model estimates
only a single set of influence parameters. This is appropriate for
networks such as undirected trade, co-sponsorship, or alliance
formation.

``` r
Y_sym = array(0, dim = c(m, m, T_len))
for (t in 1:T_len) {
    Y_t = matrix(rpois(m * m, 2), m, m)
    Y_sym[,,t] = (Y_t + t(Y_t)) / 2
    diag(Y_sym[,,t]) = NA
}
X_sym = array(0, dim = c(m, m, T_len))
for (t in 2:T_len) {
    X_sym[,,t] = Y_sym[,,t-1]
    X_sym[,,t][is.na(X_sym[,,t])] = 0
}

fit_sym = sir(
    Y = Y_sym, W = W, X = X_sym, Z = Z,
    family = "poisson", symmetric = TRUE, calc_se = TRUE
)
fit_sym
#>                       Estimate Std. Err
#> (Z) distance            0.0206   0.0315
#> (alphaW) proximity      0.0388   0.0056
#> (alphaW) shared_group   0.0186   0.0101
```

## Bipartite networks

When the sender and receiver sets differ (rectangular $\mathbf{Y}$),
bipartite structure is detected automatically and `fix_receiver = TRUE`
is enforced. This is the appropriate specification for networks where
the row and column actors represent distinct populations (e.g.,
countries sending aid to NGOs, or legislators co-sponsoring bills
introduced by other legislators). For square arrays where senders and
receivers are nonetheless distinct populations, set `bipartite = TRUE`
explicitly.

``` r
n1 = 10; n2 = 15
Y_bp = array(rpois(n1 * n2 * T_len, 2), dim = c(n1, n2, T_len))
W_bp = array(rnorm(n1 * n1 * p), dim = c(n1, n1, p))
X_bp = array(0, dim = c(n1, n2, T_len))
for (t in 2:T_len) X_bp[,,t] = log(Y_bp[,,t-1] + 1)
Z_bp = array(rnorm(n1 * n2 * 1 * T_len), dim = c(n1, n2, 1, T_len))

fit_bp = sir(
    Y = Y_bp, W = W_bp, X = X_bp, Z = Z_bp,
    family = "poisson", calc_se = TRUE
)
fit_bp
#>             Estimate Std. Err
#> (Z) Z1        0.0256   0.0238
#> (alphaW) W1  -0.0858   0.0065
#> (alphaW) W2  -0.0004   0.0078
```

## Dynamic influence covariates

In many empirical settings, the factors that mediate influence change
over time. Alliance structures evolve, trade relationships shift, and
geographic proximity may be modulated by time-varying factors such as
the presence of peacekeeping missions. When influence covariates change
over time, supply `W` as a 4D array ($m \times m \times p \times T$).
The influence parameters $\mathbf{α}$ and $\mathbf{β}$ are still
estimated as time-invariant quantities, but the influence matrices
$\mathbf{A}_{t}$ and $\mathbf{B}_{t}$ now vary with $t$ because the
underlying covariates do.

``` r
Y_dyn = array(rpois(m * m * T_len, 2), dim = c(m, m, T_len))
for (t in 1:T_len) diag(Y_dyn[,,t]) = NA

X_dyn = array(0, dim = c(m, m, T_len))
for (t in 2:T_len) X_dyn[,,t] = log(Y_dyn[,,t-1] + 1)
X_dyn[is.na(X_dyn)] = 0

# 4D W: influence covariates that change over time
W_dyn = array(rnorm(m * m * p * T_len), dim = c(m, m, p, T_len))

fit_dyn = sir(
    Y = Y_dyn, W = W_dyn, X = X_dyn,
    family = "poisson", fix_receiver = TRUE,
    calc_se = FALSE, max_iter = 10
)

# A is now a 3D array (m x m x T)
dim(fit_dyn$A)
#> [1] 15 15 10
fit_dyn$dynamic_W
#> [1] TRUE
```

## Data preparation

### Converting edge lists to arrays

Network data often arrives as edge lists (long-format data frames with
sender, receiver, time, and variable columns). The
[`cast_array()`](https://netify-dev.github.io/sir/reference/cast_array.md)
function converts this format to the multidimensional arrays required by
[`sir()`](https://netify-dev.github.io/sir/reference/sir.md).

``` r
edge_list = expand.grid(
    i = paste0("n", 1:5),
    j = paste0("n", 1:5),
    t = 1:3
)
edge_list = edge_list[edge_list$i != edge_list$j, ]
edge_list$conflict = rpois(nrow(edge_list), lambda = 2)

Y_from_el = cast_array(edge_list, var = "conflict")
dim(Y_from_el)
#> [1] 5 5 3
dimnames(Y_from_el)[[1]]
#> [1] "n1" "n2" "n3" "n4" "n5"
dimnames(Y_from_el)[[3]]
#> [1] "1" "2" "3"
```

### Constructing relational covariates

The
[`rel_covar()`](https://netify-dev.github.io/sir/reference/rel_covar.md)
function constructs relational covariate arrays from a base dyadic
variable, creating main ($z_{ij}$), reciprocal ($z_{ji}$), and
transitive ($\sum_{k}s_{ik}s_{kj}$ where $\mathbf{S}$ is the symmetrized
network) effects. These derived quantities capture common relational
mechanisms (reciprocity, transitivity) and can be included directly in
the exogenous covariate array $\mathbf{Z}$.

``` r
# base trade array (m x m x T)
trade = array(abs(rnorm(m * m * T_len)), dim = c(m, m, T_len))

# create all three relational effects
Z_trade = rel_covar(trade, "trade")
dim(Z_trade)
#> [1] 15 15  3 10
dimnames(Z_trade)[[3]]
#> [1] "trade"       "trade_recip" "trade_trans"
```

### Simulating SIR data

The [`sim_sir()`](https://netify-dev.github.io/sir/reference/sim_sir.md)
function generates synthetic network data from the SIR data-generating
process. This is useful for simulation studies, power analysis, and
verifying that the estimation procedure recovers known parameter values
under controlled conditions.

``` r
dat = sim_sir(m = 10, T_len = 8, p = 2, q = 1,
    family = "poisson", seed = 42)

str(dat[c("Y", "W", "alpha", "beta", "theta")])
#> List of 5
#>  $ Y    : num [1:10, 1:10, 1:8] 0 1 1 1 0 2 2 0 3 1 ...
#>  $ W    : num [1:10, 1:10, 1:2] 0.8895 0.0539 0.058 -0.0658 0.1555 ...
#>  $ alpha: num [1:2] 1 0.411
#>  $ beta : num [1:2] -0.169 0.109
#>  $ theta: num 0.237
```
