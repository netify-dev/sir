## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
	collapse  = TRUE,
	comment   = "#>",
	fig.align = "center",
	fig.width = 7, fig.height = 5,
	message   = FALSE, warning = FALSE
)
library(sir)
library(ggplot2)

## ----data---------------------------------------------------------------------
set.seed(42)
m = 15; T_len = 8; p = 2

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

## ----normal-------------------------------------------------------------------
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

## ----binomial-----------------------------------------------------------------
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

## ----symmetric----------------------------------------------------------------
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

## ----bipartite----------------------------------------------------------------
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

## ----dynamic-w----------------------------------------------------------------
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
fit_dyn$dynamic_W

## ----cast-array---------------------------------------------------------------
edge_list = expand.grid(
	i = paste0("n", 1:5),
	j = paste0("n", 1:5),
	t = 1:3
)
edge_list = edge_list[edge_list$i != edge_list$j, ]
edge_list$conflict = rpois(nrow(edge_list), lambda = 2)

Y_from_el = cast_array(edge_list, var = "conflict")
dim(Y_from_el)
dimnames(Y_from_el)[[1]]
dimnames(Y_from_el)[[3]]

## ----rel-covar----------------------------------------------------------------
# base trade array (m x m x T)
trade = array(abs(rnorm(m * m * T_len)), dim = c(m, m, T_len))

# create all three relational effects
Z_trade = rel_covar(trade, "trade")
dim(Z_trade)
dimnames(Z_trade)[[3]]

## ----sim-sir------------------------------------------------------------------
dat = sim_sir(m = 10, T_len = 8, p = 2, q = 1,
	family = "poisson", seed = 42)

str(dat[c("Y", "W", "alpha", "beta", "theta")])

