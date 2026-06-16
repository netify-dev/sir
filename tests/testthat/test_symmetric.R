# symmetric = TRUE now fits the genuine undirected operator A = B (the quadratic
# form A X A'); the legacy B = I upper-triangle proxy is reachable via
# fix_receiver = TRUE. these tests cover the new semantics + the migration path.

test_that("SIR fits the genuine symmetric (A = B) operator", {
	set.seed(42)
	m = 10
	T_len = 5
	p = 2
	q = 1

	# build symmetric count data
	Y = array(0, dim = c(m, m, T_len))
	for (t in 1:T_len) {
		Y_t = matrix(0, m, m)
		upper = upper.tri(Y_t)
		Y_t[upper] = rpois(sum(upper), 2)
		Y_t = Y_t + t(Y_t)
		diag(Y_t) = NA
		Y[,,t] = Y_t
	}

	W = array(rnorm(m * m * p), dim = c(m, m, p))
	for (k in seq_len(p)) {
		W[, , k] = (W[, , k] + t(W[, , k])) / 2
		diag(W[, , k]) = 0
	}
	X = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) X[,,t] = log(Y[,,t-1] + 1)
	X[is.na(X)] = 0

	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
	for (t in seq_len(T_len)) {
		for (k in seq_len(q)) {
			Z[, , k, t] = (Z[, , k, t] + t(Z[, , k, t])) / 2
		}
	}

	model = suppressWarnings(sir(Y = Y, W = W, X = X, Z = Z,
				 family = "poisson", symmetric = TRUE, seed = 1))

	expect_s3_class(model, "sir")
	expect_true(model$symmetric)
	expect_identical(model$operator, "symmetric")
	# the new operator is two-sided, NOT the legacy B = I proxy
	expect_false(model$fix_receiver)
	# A and B are the same shared operator
	expect_equal(model$A, model$B, tolerance = 1e-12)
	# diag(A) must be zero (no self-influence in the quadratic form)
	expect_lt(max(abs(diag(model$A))), 1e-10)
	# coefficient labels carry the gamma tag
	expect_true(any(grepl("\\(gammaW\\)", rownames(model$summ))))

	# fitted values are exactly symmetric
	fv = fitted(model)
	for (t in 1:T_len) {
		fv_t = fv[,,t]
		fv_t[is.na(fv_t)] = 0
		expect_equal(fv_t, t(fv_t), tolerance = 1e-10)
	}

	# predict in-sample equals fitted (shared eta helper)
	pin = predict(model)
	expect_equal(as.numeric(pin), as.numeric(model$fitted.values), tolerance = 1e-10)

	# newdata predictions stay symmetric with NA diagonal
	pred_link = predict(model, newdata = list(W = W, X = X, Z = Z), type = "link")
	for (t in 1:T_len) {
		link_t = pred_link[, , t]
		link_t[is.na(link_t)] = 0
		expect_equal(link_t, t(link_t), tolerance = 1e-10)
		expect_true(all(is.na(diag(pred_link[, , t]))))
	}

	expect_no_error(capture.output(print(model)))
	s = summary(model)
	expect_true(s$symmetric)
})

test_that("symmetric recovers a known gamma and has calibrated classical SEs", {
	dat = sim_sir(m = 14, T_len = 70, p = 3, q = 1, family = "normal",
				  symmetric = TRUE, seed = 7)
	fit = suppressWarnings(sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z,
				  family = "normal", symmetric = TRUE, seed = 1))
	# gamma recovered (alpha[1] = 1 anchor, free gamma_2..p)
	expect_equal(fit$gamma, dat$alpha, tolerance = 0.1, ignore_attr = TRUE)
	# classical bread exists; default confint is the cluster-robust interval
	expect_false(is.null(fit$vcov))
	ci_default = confint(fit)
	ci_cluster = confint(fit, se.type = "cluster")
	expect_equal(ci_default, ci_cluster)
	# classical still available on request and finite
	expect_true(all(is.finite(confint(fit, se.type = "classical"))))
})

test_that("symmetric averages asymmetric normal outcomes", {
	set.seed(42)
	m = 10
	T_len = 5

	Y = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	X = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) X[,,t] = Y[,,t-1]
	W = array(rnorm(m * m * 2), dim = c(m, m, 2))
	for (k in seq_len(2)) W[, , k] = (W[, , k] + t(W[, , k])) / 2

	expect_message(
		suppressWarnings(sir(Y = Y, W = W, X = X, family = "normal",
			symmetric = TRUE, calc_se = FALSE, seed = 1)),
		"Symmetrizing"
	)
})

test_that("symmetric rejects asymmetric discrete outcomes", {
	set.seed(43)
	m = 8
	T_len = 4

	Y = array(rpois(m * m * T_len, 2), dim = c(m, m, T_len))
	X = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) X[, , t] = Y[, , t - 1]
	W = array(rnorm(m * m * 2), dim = c(m, m, 2))
	for (k in seq_len(2)) W[, , k] = (W[, , k] + t(W[, , k])) / 2

	expect_error(
		sir(Y = Y, W = W, X = X, family = "poisson", symmetric = TRUE, seed = 1),
		"must already be symmetric"
	)
})

test_that("symmetric rejects asymmetric influence covariates", {
	set.seed(44)
	m = 8
	T_len = 4

	Y = array(0, dim = c(m, m, T_len))
	for (t in seq_len(T_len)) {
		Y_t = matrix(0, m, m)
		Y_t[upper.tri(Y_t)] = rpois(sum(upper.tri(Y_t)), 2)
		Y_t = Y_t + t(Y_t)
		diag(Y_t) = NA
		Y[, , t] = Y_t
	}
	X = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) X[, , t] = log(Y[, , t - 1] + 1)
	X[is.na(X)] = 0
	# asymmetric W slice -> should be rejected
	W = array(rnorm(m * m * 2), dim = c(m, m, 2))

	expect_error(
		sir(Y = Y, W = W, X = X, family = "poisson", symmetric = TRUE, seed = 1),
		"must be symmetric"
	)
})

test_that("symmetric is incompatible with bipartite (non-square) Y", {
	Y = array(1, dim = c(5, 6, 3))
	expect_error(
		sir(Y = Y, family = "poisson", symmetric = TRUE),
		"square"
	)
})

test_that("legacy B = I upper-triangle proxy still works via fix_receiver", {
	set.seed(42)
	m = 10; T_len = 5; p = 2
	Y = array(0, dim = c(m, m, T_len))
	for (t in 1:T_len) {
		Y_t = matrix(0, m, m)
		Y_t[upper.tri(Y_t)] = rpois(sum(upper.tri(Y_t)), 2)
		Y_t = Y_t + t(Y_t)
		diag(Y_t) = NA
		Y[, , t] = Y_t
	}
	W = array(rnorm(m * m * p), dim = c(m, m, p))
	for (k in seq_len(p)) { W[, , k] = (W[, , k] + t(W[, , k])) / 2; diag(W[, , k]) = 0 }
	X = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) X[, , t] = log(Y[, , t - 1] + 1)
	X[is.na(X)] = 0

	# the old "symmetric" behaviour (B = I, sender-side) is now fix_receiver
	model = sir(Y = Y, W = W, X = X, family = "poisson",
				fix_receiver = TRUE, calc_se = FALSE, seed = 1)
	expect_true(model$fix_receiver)
	expect_false(isTRUE(model$symmetric))
	expect_s3_class(model, "sir")
})
