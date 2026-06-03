## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
	collapse  = TRUE,
	comment   = "#>",
	fig.align = "center",
	fig.width = 7, fig.height = 5,
	message   = FALSE, warning = TRUE
)
library(sir)


## ----data---------------------------------------------------------------------
set.seed(7)
dat <- sim_sir(m = 14, T_len = 80, p = 2, q = 2, family = "poisson", seed = 7)
fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z,
           family = "poisson", calc_se = TRUE, seed = 1)


## ----classical----------------------------------------------------------------
sqrt(diag(vcov(fit, type = "classical")))   # the SEs reported by summary()
confint(fit, se.type = "classical")         # Hessian-only Wald intervals


## ----robust-------------------------------------------------------------------
se_classical <- sqrt(diag(vcov(fit, type = "classical")))
se_robust    <- sqrt(diag(vcov(fit, type = "robust")))

data.frame(
	term      = names(coef(fit)),
	classical = round(unname(se_classical), 4),
	robust    = round(unname(se_robust), 4),
	ratio     = round(unname(se_robust / se_classical), 2),
	row.names = NULL
)

confint(fit, se.type = "robust")


## ----cluster------------------------------------------------------------------
se_cluster <- sqrt(diag(vcov(fit)))

data.frame(
	term      = names(coef(fit)),
	classical = round(unname(se_classical), 4),
	cluster   = round(unname(se_cluster), 4),
	ratio     = round(unname(se_cluster / se_classical), 2),
	row.names = NULL
)


## ----cluster-confint----------------------------------------------------------
confint(fit)


## ----boot---------------------------------------------------------------------
br <- boot_sir(fit, R = 20, type = "block", seed = 123, trace = FALSE)
br_dyad <- boot_sir(fit, type = "dyad", seed = 123, trace = FALSE)

data.frame(
	method = c("Block Bootstrap", "Dyad Jackknife"),
	valid_refits = c(br$n_valid, br_dyad$n_valid),
	total_refits = c(br$n_total, br_dyad$n_total),
	interval_type = c(br$interval, br_dyad$interval)
)


## ----confint-block------------------------------------------------------------
confint(fit, boot = br)


## ----confint-dyad-------------------------------------------------------------
confint(fit, boot = br_dyad)


## ----reporting-table----------------------------------------------------------
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


## ----fix-receiver-------------------------------------------------------------
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


## ----compare------------------------------------------------------------------
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


## ----forecast-cv--------------------------------------------------------------
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

