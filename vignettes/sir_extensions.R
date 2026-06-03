## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
	collapse  = TRUE,
	comment   = "#>",
	fig.align = "center",
	fig.width = 7, fig.height = 5,
	message   = FALSE, warning = TRUE
)
library(sir)


## ----normal-------------------------------------------------------------------
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


## ----normal-response----------------------------------------------------------
mu_norm <- predict(fit_norm)
round(quantile(mu_norm[!is.na(mu_norm)], c(0.05, 0.5, 0.95)), 3)


## ----binomial-----------------------------------------------------------------
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


## ----binomial-prob------------------------------------------------------------
mu_bin <- predict(fit_bin)                 # fitted tie probabilities
round(quantile(mu_bin[!is.na(mu_bin)], c(0.05, 0.5, 0.95)), 3)


## ----symmetric----------------------------------------------------------------
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


## ----bipartite----------------------------------------------------------------
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


## ----bipartite-full-----------------------------------------------------------
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

data.frame(
	term = names(coef(fit_full)),
	estimate = round(unname(coef(fit_full)), 3),
	target = c(0.6, 0.8, -0.5),
	row.names = NULL
)


## ----bipartite-full-boot------------------------------------------------------
bj <- boot_sir(fit_full, type = "dyad", seed = 1)
confint(bj)


## ----dynamic-w----------------------------------------------------------------
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
fit_dyn$dynamic_W


## ----dynamic-w-edges----------------------------------------------------------
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


## ----cast-array---------------------------------------------------------------
set.seed(4)
edge_list <- expand.grid(i = paste0("n", 1:5), j = paste0("n", 1:5), t = 1:3)
edge_list <- edge_list[edge_list$i != edge_list$j, ]
edge_list$conflict <- rpois(nrow(edge_list), lambda = 2)

Y_from_el <- cast_array(edge_list, var = "conflict")
dim(Y_from_el)
dimnames(Y_from_el)[[1]]


## ----rel-covar----------------------------------------------------------------
set.seed(5)
trade   <- array(abs(rnorm(8 * 8 * 4)), dim = c(8, 8, 4))   # base dyadic variable
Z_trade <- rel_covar(trade, "trade")
dim(Z_trade)
dimnames(Z_trade)[[3]]


## ----sim-sir------------------------------------------------------------------
dat <- sim_sir(m = 10, T_len = 8, p = 2, q = 1, family = "poisson", seed = 42)
str(dat[c("Y", "W", "X", "Z", "alpha", "beta", "theta")], max.level = 1)

