## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
	collapse  = TRUE,
	comment   = "#>",
	fig.align = "center",
	fig.width = 7, fig.height = 5,
	message   = FALSE, warning = FALSE
)
library(sir)

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

## ----fit----------------------------------------------------------------------
fit_als = sir(
	Y = Y, W = W, X = X, Z = Z,
	family = "poisson", method = "ALS", calc_se = TRUE
)

## ----vcov---------------------------------------------------------------------
V_classical = vcov(fit_als, type = "classical")
V_robust = vcov(fit_als, type = "robust")

# classical and robust standard errors
if (!is.null(V_classical)) sqrt(diag(V_classical))
if (!is.null(V_robust)) sqrt(diag(V_robust))

## ----fix-receiver-------------------------------------------------------------
fit_fix = sir(
	Y = Y, W = W, X = X, Z = Z,
	family = "poisson",
	method = "ALS",
	fix_receiver = TRUE,
	calc_se = TRUE
)

summary(fit_fix)

## ----compare------------------------------------------------------------------
AIC(fit_als)
AIC(fit_fix)
BIC(fit_als)
BIC(fit_fix)

## ----boot, cache = TRUE-------------------------------------------------------
boot_results = boot_sir(
	fit_als,
	R = 100,
	type = "block",
	seed = 123,
	trace = FALSE
)

boot_results

## ----confint------------------------------------------------------------------
confint(fit_als)

## ----confint-boot-------------------------------------------------------------
confint(fit_als, boot = boot_results)

