# check cluster-robust variance calculations

test_that("per-observation scores reconstruct the kernel score outer-product", {
	set.seed(1)
	dat = sim_sir(m = 10, T_len = 20, p = 2, q = 2, family = "poisson", seed = 1)
	fit = sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z,
			  family = "poisson", calc_se = TRUE, seed = 1)

	sc = sir:::.sir_obs_scores(fit)
	Z_list = sir:::prepare_Z_list(dat$Z)
	gH = sir:::cpp_mll_gH(fit$tab, dat$Y, dat$W, dat$X, Z_list, "poisson")

	# crossprod(scores) == kernel shess (the BHHH meat); colSums(scores) == -grad.
	expect_equal(dim(crossprod(sc$scores)), dim(gH$shess))
	expect_lt(max(abs(crossprod(sc$scores) - gH$shess)) / max(abs(gH$shess)), 1e-8)
	expect_lt(max(abs(colSums(sc$scores) + as.numeric(gH$grad))), 1e-5)
})

test_that("cluster-robust vcov is a valid PSD covariance, exposed via vcov/confint", {
	set.seed(2)
	dat = sim_sir(m = 10, T_len = 20, p = 2, q = 2, family = "poisson", seed = 2)
	fit = sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z,
			  family = "poisson", calc_se = TRUE, seed = 2)

	V = vcov(fit, type = "cluster")
	expect_equal(dim(V), dim(fit$vcov))
	expect_true(all(is.finite(V)))
	expect_gt(min(eigen(V, symmetric = TRUE, only.values = TRUE)$values), -1e-8)
	expect_true(isSymmetric(unname(V), tol = 1e-8))

	ci = confint(fit, se.type = "cluster")
	expect_equal(nrow(ci), length(coef(fit)))
	expect_true(all(is.finite(ci)))
})

# dynamic (4D) W: cluster-robust SEs are first-class. the per-cell scores match
# the dynamic c++ kernel's gradient at the fitted tab (the score equation), and
# the sandwich is a valid PSD covariance that the default accessors route to.

make_dyn_W = function(m, T_len, p, seed, sym = FALSE) {
	set.seed(seed)
	W = array(0, c(m, m, p, T_len))
	W[, , , 1] = rnorm(m * m * p)
	for (t in 2:T_len) {
		W[, , , t] = 0.9 * W[, , , t - 1] +
			sqrt(1 - 0.81) * array(rnorm(m * m * p), c(m, m, p))
	}
	for (t in seq_len(T_len)) for (k in seq_len(p)) {
		if (sym) { a = W[, , k, t]; W[, , k, t] = (a + t(a)) / 2 }
		diag(W[, , k, t]) = 0
	}
	W
}

test_that("full-directed dynamic scores match the dynamic c++ kernel gradient", {
	m = 10; T_len = 35; p = 2
	W = make_dyn_W(m, T_len, p, seed = 31)
	A_dyn = function(t) W[, , 1, t] + 0.5 * W[, , 2, t]
	B_dyn = function(t) W[, , 1, t] - 0.3 * W[, , 2, t]
	set.seed(5)
	Z = array(rnorm(m * m * T_len), c(m, m, T_len))
	X = array(rnorm(m * m * T_len), c(m, m, T_len)) / (m - 1)
	Y = array(0, c(m, m, T_len))
	for (t in seq_len(T_len)) {
		Y[, , t] = 0.5 * Z[, , t] + A_dyn(t) %*% X[, , t] %*% t(B_dyn(t)) +
			matrix(rnorm(m * m, sd = 0.3), m, m)
		diag(Y[, , t]) = NA
	}
	fit = suppressWarnings(sir(Y, W = W, X = X, Z = Z, family = "normal", seed = 1))
	expect_null(fit$se_source)

	# the per-cell scores reproduce the dynamic kernel's gradient at the fitted
	# tab: colSums(scores) == -grad (the score equation). use the fit's own X/W/Z,
	# which sir() stores with the self-tie diagonal zeroed.
	sc = sir:::.sir_obs_scores(fit)
	W_field = sir:::prepare_W_field(fit$W)
	Z_list = sir:::prepare_Z_list(fit$Z)
	gH = sir:::cpp_mll_gH_dyn(fit$tab, fit$Y, W_field, fit$X, Z_list, "normal")
	expect_lt(max(abs(colSums(sc$scores) + as.numeric(gH$grad))), 1e-3)

	# the sandwich is PSD with a t(G - 1) reference and is the default
	V = vcov(fit, type = "cluster")
	expect_gt(min(eigen(V, symmetric = TRUE, only.values = TRUE)$values), -1e-8)
	expect_equal(attr(V, "cluster_df"), m - 1L)
	expect_equal(vcov(fit), V)
	se_classical = sqrt(diag(vcov(fit, type = "classical")))
	expect_false(isTRUE(all.equal(se_classical, sqrt(diag(V)))))
})

test_that("dynamic (4D) W gets a valid cluster-robust sandwich by default", {
	m = 10; T_len = 45; p = 2
	W = make_dyn_W(m, T_len, p, seed = 11)
	A_dyn = function(t) W[, , 1, t] + 0.7 * W[, , 2, t]
	Y = array(0, c(m, m, T_len)); X = array(0, c(m, m, T_len))
	for (t in seq_len(T_len)) {
		if (t > 1) { lg = Y[, , t - 1]; lg[is.na(lg)] = 0; X[, , t] = lg / (m - 1) }
		Y[, , t] = A_dyn(t) %*% X[, , t] + matrix(rnorm(m * m, sd = 0.4), m, m)
		diag(Y[, , t]) = NA
	}
	fit = suppressWarnings(sir(Y, W = W, X = X, family = "normal",
							   fix_receiver = TRUE, seed = 1))
	# analytic SEs are available, so the jackknife fallback did not fire
	expect_null(fit$se_source)

	# the reconstructed per-cell scores satisfy the MLE score equation
	sc = sir:::.sir_obs_scores(fit)
	expect_lt(max(abs(colSums(sc$scores))), 1e-4)

	# the sandwich is a valid PSD covariance with a t(G - 1) reference (G = m)
	V = vcov(fit, type = "cluster")
	expect_equal(dim(V), dim(fit$vcov))
	expect_true(all(is.finite(V)))
	expect_gt(min(eigen(V, symmetric = TRUE, only.values = TRUE)$values), -1e-8)
	expect_equal(attr(V, "cluster_df"), m - 1L)

	# the default vcov()/confint()/tidy() now route to cluster, not classical
	expect_equal(vcov(fit), V)
	expect_true(all(is.finite(confint(fit))))
	expect_true(all(is.finite(tidy(fit)$std.error)))

	# the sandwich actually differs from the classical SEs it replaces
	se_classical = sqrt(diag(vcov(fit, type = "classical")))
	se_cluster = sqrt(diag(V))
	expect_false(isTRUE(all.equal(se_classical, se_cluster)))
})

test_that("symmetric dynamic (4D) W gets a cluster-robust sandwich", {
	m = 10; T_len = 45; p = 2
	W = make_dyn_W(m, T_len, p, seed = 21, sym = TRUE)
	A_dyn = function(t) 0.9 * W[, , 1, t] + 0.5 * W[, , 2, t]
	X = array(0, c(m, m, T_len)); Y = array(0, c(m, m, T_len))
	for (t in seq_len(T_len)) {
		a = matrix(rnorm(m * m), m, m); X[, , t] = (a + t(a)) / 2 / (m - 1)
		e = A_dyn(t) %*% X[, , t] %*% t(A_dyn(t))
		y = e + matrix(rnorm(m * m, sd = 0.4), m, m); y = (y + t(y)) / 2
		diag(y) = NA; Y[, , t] = y
	}
	fit = suppressWarnings(sir(Y, W = W, X = X, family = "normal",
							   symmetric = TRUE, seed = 1))
	expect_null(fit$se_source)

	# symmetric ALS converges block-coordinate, so check the score equation
	# relative to the score magnitude rather than to machine zero
	sc = sir:::.sir_obs_scores(fit)
	expect_lt(max(abs(colSums(sc$scores))) / max(abs(sc$scores)), 1e-2)

	V = vcov(fit, type = "cluster")
	expect_true(all(is.finite(V)))
	expect_equal(attr(V, "cluster_df"), m - 1L)
	expect_true(all(is.finite(confint(fit))))
})
