# m4: optim-vs-ALS agreement (fix_receiver to avoid bilinear local optima)
# m5: sign-consistency across seeds (fix_receiver, stronger signal)
# m10: predict round-trip equals fitted

test_that("ALS and optim produce similar estimates (fix_receiver)", {
	set.seed(42)
	m = 10
	T_len = 8
	p = 2
	q = 1

	Y = array(rpois(m * m * T_len, 2), dim = c(m, m, T_len))
	for (t in 1:T_len) diag(Y[,,t]) = NA
	W = array(rnorm(m * m * p, sd = 0.2), dim = c(m, m, p))
	X = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) { X[,,t] = Y[,,t-1]; X[is.na(X[,,t])] = 0 }
	Z = array(rnorm(m * m * q * T_len, sd = 0.3), dim = c(m, m, q, T_len))

	fit_als = sir(Y = Y, W = W, X = X, Z = Z, family = "poisson",
				method = "ALS", fix_receiver = TRUE,
				calc_se = FALSE, trace = FALSE, max_iter = 100)
	fit_opt = sir(Y = Y, W = W, X = X, Z = Z, family = "poisson",
				method = "optim", fix_receiver = TRUE,
				calc_se = FALSE, trace = 0)

	# both should produce finite estimates
	expect_true(all(is.finite(fit_als$tab)))
	expect_true(all(is.finite(fit_opt$tab)))

	# log-likelihoods should be close (within 5%)
	ll_als = as.numeric(logLik(fit_als))
	ll_opt = as.numeric(logLik(fit_opt))
	expect_true(abs(ll_als - ll_opt) / abs(ll_als) < 0.05)

	# parameter estimates should be similar
	expect_equal(unname(fit_als$tab), unname(fit_opt$tab), tolerance = 0.5)
})

test_that("sign of influence parameters is consistent across seeds", {
	m = 12
	T_len = 15
	p = 2

	# fix_receiver model: tab = [alpha_1, alpha_2]
	true_alpha = c(0.5, 0.3)

	set.seed(100)
	W = array(rnorm(m * m * p, sd = 0.2), dim = c(m, m, p))

	signs = rep(NA, 5)

	for (s in 1:5) {
	set.seed(s * 100)
	X = array(0, dim = c(m, m, T_len))
	Y = array(0, dim = c(m, m, T_len))

	for (t in 1:T_len) {
		if (t == 1) {
		X[,,t] = matrix(rpois(m * m, 2), m, m)
		diag(X[,,t]) = 0
		} else {
		X[,,t] = Y[,,t-1]
		X[is.na(X[,,t])] = 0
		}
		ETA_t = eta_tab(true_alpha, W, X[,,t, drop=FALSE], Z = NULL, fix_receiver = TRUE)
		lambda_t = exp(ETA_t[,,1])
		diag(lambda_t) = 0
		lambda_t = pmin(pmax(lambda_t, 0.01), 50)
		Y_t = matrix(rpois(m * m, lambda = c(lambda_t)), m, m)
		diag(Y_t) = NA
		Y[,,t] = Y_t
	}

	X_fit = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) { X_fit[,,t] = Y[,,t-1]; X_fit[is.na(X_fit[,,t])] = 0 }

	model = sir(Y = Y, W = W, X = X_fit,
				family = "poisson", fix_receiver = TRUE,
				calc_se = FALSE, trace = FALSE, max_iter = 50)
	est = unname(model$tab)
	signs[s] = sign(est[2])  # alpha_2
	}

	# at least 4 out of 5 seeds should agree on the sign
	expect_true(max(table(signs)) >= 4, label = "alpha_2 sign consistency")
})

test_that("predict(model) equals fitted(model) for training data", {
	set.seed(42)
	m = 8
	T_len = 3
	p = 2
	q = 1

	Y = array(rpois(m * m * T_len, 2), dim = c(m, m, T_len))
	for (t in 1:T_len) diag(Y[,,t]) = NA
	W = array(rnorm(m * m * p, sd = 0.2), dim = c(m, m, p))
	X = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) { X[,,t] = Y[,,t-1]; X[is.na(X[,,t])] = 0 }
	Z = array(rnorm(m * m * q * T_len, sd = 0.3), dim = c(m, m, q, T_len))

	model = sir(Y = Y, W = W, X = X, Z = Z, family = "poisson",
				calc_se = FALSE, trace = FALSE, max_iter = 20)

	pred = predict(model, type = "response")
	fit = fitted(model)

	expect_equal(pred, fit)
})
