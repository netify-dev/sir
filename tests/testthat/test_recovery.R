
# recovery tests for the influence structure.
# the identified object is the rank-at-most-one matrix C = alpha %*% t(beta), invariant
# to the reciprocal-scale nuisance (c*alpha, beta/c), so recovery is measured by
# the relative frobenius error relC = ||C_hat - C_true||_F / ||C_true||_F, taken
# as the median over several seeds for a stable measure.

relC_of_C = function(fit, alpha_true, beta_true) {
	p = length(alpha_true)
	q = length(fit$theta)
	alpha_hat = c(1, if (p > 1) fit$tab[(q + 1):(q + p - 1)] else numeric(0))
	beta_hat  = fit$tab[(q + p):(q + 2 * p - 1)]
	C_hat = outer(alpha_hat, beta_hat)
	C_true = outer(alpha_true, beta_true)
	sqrt(sum((C_hat - C_true)^2)) / sqrt(sum(C_true^2))
}

# median relC over `nseed` independent simulated datasets (sim_sir's default
# modest, stationary coefficients).
median_relC = function(family, T_len, p = 2, nseed = 5, base = 100,
						m = 14, alpha = NULL, beta = NULL) {
	rc = vapply(seq_len(nseed), function(s) {
		dat = sim_sir(m = m, T_len = T_len, p = p, q = 2, family = family,
					   alpha = alpha, beta = beta, seed = base + s)
		fit = tryCatch(
			sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = family,
				method = "ALS", calc_se = FALSE, seed = 1),
			error = function(e) NULL)
		if (is.null(fit)) return(NA_real_)
		relC_of_C(fit, dat$alpha, dat$beta)
	}, numeric(1))
	median(rc, na.rm = TRUE)
}

test_that("Poisson influence coefficients (C = alpha beta') recover", {
	expect_lt(median_relC("poisson", T_len = 80), 0.40)
})

test_that("Normal influence coefficients recover", {
	expect_lt(median_relC("normal", T_len = 80), 0.40)
})

test_that("Binomial influence coefficients recover (with identifiable signal)", {
	# binary outcomes have a high noise floor, so recovery needs a larger influence
	# signal: stronger coefficients and a slightly larger network
	expect_lt(
		median_relC("binomial", T_len = 300, m = 20,
					alpha = c(1, 1.6), beta = c(1.4, -1.3)),
		0.45)
})

test_that("influence recovery holds with p = 3 influence covariates", {
	expect_lt(median_relC("poisson", T_len = 80, p = 3), 0.45)
})

test_that("direct effects (theta) recover tightly", {
	errs = vapply(1:5, function(s) {
		dat = sim_sir(m = 14, T_len = 80, p = 2, q = 2, family = "poisson", seed = 200 + s)
		fit = sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson",
				   calc_se = FALSE, seed = 1)
		max(abs(fit$theta - dat$theta))
	}, numeric(1))
	expect_lt(median(errs), 0.2)
})

test_that("recovery is far below the broken near-identity-anchor regime", {
	expect_lt(median_relC("poisson", T_len = 120), 0.45)
})

# the simulator's linear predictor must equal eta_tab() at the true parameters
test_that("sim_sir eta matches eta_tab() at the true parameters", {
	set.seed(42)
	dat = sim_sir(m = 12, T_len = 40, p = 2, q = 2, family = "poisson", seed = 42)
	tab_true = c(dat$theta, dat$alpha[-1], dat$beta)
	eta_est = sir:::eta_tab(tab_true, dat$W, dat$X, dat$Z)

	eta_sim = array(0, dim = dim(dat$X))
	for (t in seq_len(dim(dat$X)[3])) eta_sim[, , t] = dat$A %*% dat$X[, , t] %*% t(dat$B)
	for (k in seq_len(dim(dat$Z)[3])) {
		for (t in seq_len(dim(dat$X)[3])) {
			eta_sim[, , t] = eta_sim[, , t] + dat$theta[k] * dat$Z[, , k, t]
		}
	}

	# compare off-diagonal only (the diagonal/self-tie is not modeled)
	mask = array(TRUE, dim = dim(eta_sim))
	for (t in seq_len(dim(mask)[3])) diag(mask[, , t]) = FALSE
	expect_lt(max(abs(eta_est - eta_sim)[mask]), 1e-8)
})

# the unbounded-family recursion must stay stationary across seeds
test_that("sim_sir Poisson stays stationary across seeds (no count explosion)", {
	ymax = vapply(201:215, function(s) {
		dat = sim_sir(m = 14, T_len = 80, p = 2, q = 2, family = "poisson", seed = s)
		max(dat$Y)
	}, numeric(1))
	expect_true(all(ymax < 1e5))
})
