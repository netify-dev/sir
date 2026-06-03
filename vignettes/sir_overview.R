## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
	collapse  = TRUE,
	comment   = "#>",
	fig.align = "center",
	fig.width = 7, fig.height = 5,
	message   = FALSE, warning = TRUE
)
library(sir)
qval <- function(x, prob) unname(stats::quantile(x, prob, na.rm = TRUE))


## ----simulate-----------------------------------------------------------------
set.seed(42)
dat <- sim_sir(m = 14, T_len = 80, p = 2, q = 2, family = "poisson", seed = 42)

fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 1)
fit


## ----recovery-----------------------------------------------------------------
truth     <- c(dat$theta, dat$alpha[-1], dat$beta)
estimated <- coef(fit)

data.frame(
	term      = names(estimated),
	truth     = round(truth, 3),
	estimate  = round(unname(estimated), 3),
	abs_error = round(abs(truth - unname(estimated)), 3)
)


## ----recovery-c---------------------------------------------------------------
C_truth <- outer(dat$alpha, dat$beta)
C_estimated <- outer(fit$alpha, fit$beta)

data.frame(
	quantity = c("maximum absolute C error", "relative Frobenius C error"),
	value = round(c(
		max(abs(C_truth - C_estimated)),
		sqrt(sum((C_truth - C_estimated)^2)) / sqrt(sum(C_truth^2))
	), 3)
)


## ----coefficient-table--------------------------------------------------------
ci <- confint(fit)
data.frame(
	term = names(coef(fit)),
	estimate = round(unname(coef(fit)), 3),
	cluster_low = round(ci[, 1], 3),
	cluster_high = round(ci[, 2], 3),
	row.names = NULL
)


## ----influence-cells----------------------------------------------------------
A <- fit$A
diag(A) <- NA                                      # self-influence is not modeled
ord <- order(abs(A), decreasing = TRUE, na.last = NA)[1:5]
data.frame(
	source_node     = ((ord - 1) %/% nrow(A)) + 1,   # column k (the influencer)
	influenced_node = ((ord - 1) %%  nrow(A)) + 1,   # row i (the influenced)
	influence       = round(A[ord], 3)
)


## ----receiver-cells-----------------------------------------------------------
B <- fit$B
diag(B) <- NA
ord_b <- order(abs(B), decreasing = TRUE, na.last = NA)[1:5]
data.frame(
	source_receiver     = ((ord_b - 1) %/% nrow(B)) + 1,
	influenced_receiver = ((ord_b - 1) %%  nrow(B)) + 1,
	influence           = round(B[ord_b], 3)
)


## ----influence-magnitude------------------------------------------------------
round(quantile(A[!is.na(A)], c(0.05, 0.25, 0.5, 0.75, 0.95)), 3)


## ----plots, fig.width = 8, fig.height = 8, fig.cap = "SIR diagnostics for the simulated fit: sender and receiver influence heatmaps, followed by off-diagonal influence-weight distributions.", fig.alt = "Four diagnostic panels showing the A and B influence matrix heatmaps and the distributions of their off-diagonal entries."----
plot(fit, which = 1:4)


## ----network, fig.height = 5, fig.cap = "Network view of the strongest fitted sender-side influence channels when igraph and ggraph are installed; otherwise the chunk prints the corresponding top-channel table.", fig.alt = "Directed network plot of sender-side influence channels, with arrows from source nodes to influenced nodes."----
if (requireNamespace("ggraph", quietly = TRUE) &&
	requireNamespace("igraph", quietly = TRUE)) {
	plot_sir_network(fit, matrix = "A", threshold = 0.15)
} else {
	A_network <- fit$A
	diag(A_network) <- NA
	edge_order <- order(abs(A_network), decreasing = TRUE, na.last = NA)[1:5]
	data.frame(
		source_node = ((edge_order - 1) %/% nrow(A_network)) + 1,
		influenced_node = ((edge_order - 1) %% nrow(A_network)) + 1,
		influence = round(A_network[edge_order], 3)
	)
}


## ----predict------------------------------------------------------------------
mu_hat <- predict(fit)                             # in-sample expected counts
dim(mu_hat)

Zcf <- dat$Z
Zcf[, , 1, ] <- Zcf[, , 1, ] + sd(Zcf[, , 1, ])    # shift first Z covariate up 1 sd
mu_cf <- predict(fit, newdata = list(W = dat$W, X = dat$X, Z = Zcf))
delta_z <- mu_cf - mu_hat

data.frame(
	scenario = "Increase Z1 by 1 SD",
	baseline_mean = mean(mu_hat, na.rm = TRUE),
	scenario_mean = mean(mu_cf, na.rm = TRUE),
	mean_change = mean(delta_z, na.rm = TRUE),
	median_change = stats::median(delta_z, na.rm = TRUE),
	q05_change = qval(delta_z, 0.05),
	q95_change = qval(delta_z, 0.95),
	check.names = FALSE,
	row.names = NULL
)


## ----scenario-helper----------------------------------------------------------
W_base <- dat$W
off <- row(W_base[, , 1]) != col(W_base[, , 1])
w_step <- stats::sd(W_base[, , 1][off], na.rm = TRUE)

w_scenarios <- lapply(c("-1 SD", "baseline", "+1 SD"), function(label) {
	W_tmp <- W_base
	W1 <- W_tmp[, , 1]
	shift <- switch(label, "-1 SD" = -w_step, "baseline" = 0, "+1 SD" = w_step)
	W1[off] <- W1[off] + shift
	W_tmp[, , 1] <- W1
	mu_tmp <- predict(fit, newdata = list(W = W_tmp, X = dat$X, Z = dat$Z))
	delta <- mu_tmp - mu_hat
	data.frame(
		scenario = label,
		scenario_mean = mean(mu_tmp, na.rm = TRUE),
		mean_change = mean(delta, na.rm = TRUE),
		median_change = stats::median(delta, na.rm = TRUE),
		p90_abs_change = qval(abs(delta), 0.9),
		check.names = FALSE,
		row.names = NULL
	)
})
do.call(rbind, w_scenarios)


## ----icews--------------------------------------------------------------------
data(icews)
icews_nodes <- 1:10
icews_periods <- 1:24

Y_icews <- icews$Y[icews_nodes, icews_nodes, icews_periods, drop = FALSE]
X_icews <- icews$X[icews_nodes, icews_nodes, icews_periods, drop = FALSE]
W_icews <- icews$W[icews_nodes, icews_nodes, , drop = FALSE]
Z_icews <- icews$Z[icews_nodes, icews_nodes, , icews_periods, drop = FALSE]

data.frame(
	Component = c("Y", "X", "W", "Z"),
	Dimensions = c(
		paste(dim(Y_icews), collapse = " x "),
		paste(dim(X_icews), collapse = " x "),
		paste(dim(W_icews), collapse = " x "),
		paste(dim(Z_icews), collapse = " x ")
	)
)

ifit <- sir(
	Y_icews,
	W = W_icews,
	X = X_icews,
	Z = Z_icews,
	family = "poisson",
	seed = 1,
	max_iter = 20
)


## ----icews-diagnostics--------------------------------------------------------
data.frame(
	Diagnostic = c("ALS converged", "Classical SEs reliable"),
	Value = c(ifit$convergence, ifit$se_reliable)
)


## ----icews-coef---------------------------------------------------------------
term_labels <- c(
	"(Z) mConf" = "Direct: Lagged Material Conflict",
	"(Z) mConf_ji" = "Direct: Reciprocal Lagged Material Conflict",
	"(Z) minDistLog" = "Direct: Minimum Logged Distance",
	"(Z) ally" = "Direct: Alliance",
	"(Z) verbCoop" = "Direct: Verbal Cooperation",
	"(alphaW) ally" = "Sender Channel: Alliance",
	"(alphaW) verbCoop" = "Sender Channel: Verbal Cooperation",
	"(alphaW) minDistLog" = "Sender Channel: Minimum Logged Distance",
	"(betaW) int" = "Receiver Channel: Baseline",
	"(betaW) ally" = "Receiver Channel: Alliance",
	"(betaW) verbCoop" = "Receiver Channel: Verbal Cooperation",
	"(betaW) minDistLog" = "Receiver Channel: Minimum Logged Distance"
)

se_classical_icews <- sqrt(diag(vcov(ifit, type = "classical")))
se_cluster_icews <- sqrt(diag(vcov(ifit, type = "cluster")))

data.frame(
	"Term" = unname(term_labels[names(coef(ifit))]),
	"Estimate" = round(unname(coef(ifit)), 3),
	"Classical SE" = signif(unname(se_classical_icews), 4),
	"Cluster SE" = signif(unname(se_cluster_icews), 4),
	"Cluster/Classical" = round(unname(se_cluster_icews / se_classical_icews), 2),
	check.names = FALSE
)


## ----icews-cells--------------------------------------------------------------
A <- ifit$A
dimnames(A) <- list(icews$countries[icews_nodes], icews$countries[icews_nodes])
diag(A) <- NA
top_edges <- order(abs(A), decreasing = TRUE, na.last = NA)[1:5]
source_index <- ((top_edges - 1) %/% nrow(A)) + 1
influenced_index <- ((top_edges - 1) %% nrow(A)) + 1

data.frame(
	Source = colnames(A)[source_index],
	Influenced = rownames(A)[influenced_index],
	Weight = round(A[top_edges], 2)
)


## ----icews-b-cells------------------------------------------------------------
B <- ifit$B
dimnames(B) <- list(icews$countries[icews_nodes], icews$countries[icews_nodes])
diag(B) <- NA
top_b <- order(abs(B), decreasing = TRUE, na.last = NA)[1:5]
b_source_index <- ((top_b - 1) %/% nrow(B)) + 1
b_influenced_index <- ((top_b - 1) %% nrow(B)) + 1

data.frame(
	Source_receiver = colnames(B)[b_source_index],
	Influenced_receiver = rownames(B)[b_influenced_index],
	Weight = round(B[top_b], 2)
)


## ----icews-scenario-----------------------------------------------------------
mu_icews <- predict(ifit)
X_scen <- X_icews
src <- source_index[1]
x_shift <- 0.25 * stats::sd(X_icews[src, , ], na.rm = TRUE)
X_scen[src, , ] <- X_scen[src, , ] + x_shift
mu_scen <- predict(ifit, newdata = list(W = W_icews, X = X_scen, Z = Z_icews))
delta <- mu_scen - mu_icews

data.frame(
	scenario = paste("Increase lagged signal sent by", colnames(A)[src], "by 0.25 SD"),
	baseline_mean = mean(mu_icews, na.rm = TRUE),
	scenario_mean = mean(mu_scen, na.rm = TRUE),
	mean_change = mean(delta, na.rm = TRUE),
	median_change = stats::median(delta, na.rm = TRUE),
	q05_change = qval(delta, 0.05),
	q95_change = qval(delta, 0.95),
	p90_abs_change = qval(abs(delta), 0.9),
	max_abs_change = max(abs(delta), na.rm = TRUE),
	check.names = FALSE,
	row.names = NULL
)


## ----icews-scenario-dyads-----------------------------------------------------
top_delta <- order(abs(delta), decreasing = TRUE, na.last = NA)[1:5]
top_idx <- arrayInd(top_delta, dim(delta))
data.frame(
	Source = icews$countries[icews_nodes][top_idx[, 1]],
	Receiver = icews$countries[icews_nodes][top_idx[, 2]],
	Month = icews_periods[top_idx[, 3]],
	baseline_mean = round(mu_icews[top_delta], 3),
	scenario_mean = round(mu_scen[top_delta], 3),
	delta = round(delta[top_delta], 3),
	row.names = NULL
)

