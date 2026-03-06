# test that sir() can recover known parameters from simulated data.
# uses the model's own eta_tab to generate data, then fits and checks recovery.

test_that("sir recovers Poisson parameters from simulated data", {
	set.seed(123)
	m = 15
	T_len = 30
	p = 2
	q = 1

	# true parameters: [theta_1, alpha_2, beta_1, beta_2]
	true_theta = 0.2
	true_alpha = c(1, 0.3)   # alpha_1=1 fixed
	true_beta  = c(0.4, -0.2)
	true_tab   = c(true_theta, true_alpha[2], true_beta)

	# generate covariates (small sd to keep eta moderate)
	W = array(rnorm(m * m * p, sd = 0.15), dim = c(m, m, p))
	Z = array(rnorm(m * m * q * T_len, sd = 0.3), dim = c(m, m, q, T_len))

	X = array(0, dim = c(m, m, T_len))
	Y = array(0, dim = c(m, m, T_len))

	# simulate forward
	for (t in 1:T_len) {
	if (t == 1) {
		X[,,t] = matrix(rpois(m * m, 2), m, m)
		diag(X[,,t]) = 0
	} else {
		X[,,t] = Y[,,t-1]
		X[is.na(X[,,t])] = 0
	}
	ETA_t = sir::eta_tab(true_tab, W, X[,,t, drop=FALSE], Z[,,,t, drop=FALSE])
	lambda_t = exp(ETA_t[,,1])
	diag(lambda_t) = 0
	lambda_t = pmin(lambda_t, 50)
	lambda_t = pmax(lambda_t, 0.01)
	Y_t = matrix(rpois(m * m, lambda = c(lambda_t)), m, m)
	diag(Y_t) = NA
	Y[,,t] = Y_t
	}

	# fit using lagged Y as X
	X_fit = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) {
	X_fit[,,t] = Y[,,t-1]
	X_fit[is.na(X_fit[,,t])] = 0
	}

	model = sir(Y = Y, W = W, X = X_fit, Z = Z,
				 family = "poisson", calc_se = FALSE,
				 trace = FALSE, max_iter = 50)

	est = unname(model$tab)

	# with T=30, estimates should be close to true values
	# tab = [theta_1, alpha_2, beta_1, beta_2]
	expect_equal(est[1], true_theta, tolerance = 0.25,
				 label = "theta recovery")
	expect_equal(est[2], true_alpha[2], tolerance = 0.25,
				 label = "alpha_2 recovery")
	expect_equal(est[3], true_beta[1], tolerance = 0.25,
				 label = "beta_1 recovery")
	expect_equal(est[4], true_beta[2], tolerance = 0.25,
				 label = "beta_2 recovery")
})

test_that("sir recovers Normal parameters from simulated data", {
	set.seed(456)
	m = 15
	T_len = 8
	p = 2
	q = 1

	true_theta = -0.3
	true_alpha = c(1, 0.3)
	true_beta  = c(0.4, 0.2)
	true_tab   = c(true_theta, true_alpha[2], true_beta)

	W = array(rnorm(m * m * p, sd = 0.15), dim = c(m, m, p))
	Z = array(rnorm(m * m * q * T_len, sd = 0.3), dim = c(m, m, q, T_len))
	X = array(0, dim = c(m, m, T_len))
	Y = array(0, dim = c(m, m, T_len))

	for (t in 1:T_len) {
	if (t == 1) {
		X[,,t] = matrix(rnorm(m * m, sd = 0.5), m, m)
		diag(X[,,t]) = 0
	} else {
		X[,,t] = Y[,,t-1]
		X[is.na(X[,,t])] = 0
	}
	ETA_t = sir::eta_tab(true_tab, W, X[,,t, drop=FALSE], Z[,,,t, drop=FALSE])
	Y_t = ETA_t[,,1] + rnorm(m * m, sd = 1)
	diag(Y_t) = NA
	Y[,,t] = Y_t
	}

	X_fit = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) {
	X_fit[,,t] = Y[,,t-1]
	X_fit[is.na(X_fit[,,t])] = 0
	}

	model = sir(Y = Y, W = W, X = X_fit, Z = Z,
				 family = "normal", calc_se = FALSE,
				 trace = FALSE, max_iter = 50)

	est = unname(model$tab)

	expect_equal(est[1], true_theta, tolerance = 0.5,
				 label = "Normal theta recovery")
	expect_equal(est[2], true_alpha[2], tolerance = 0.5,
				 label = "Normal alpha_2 recovery")
	expect_equal(est[3], true_beta[1], tolerance = 0.5,
				 label = "Normal beta_1 recovery")
})

test_that("sir recovers fix_receiver parameters", {
	set.seed(789)
	m = 20
	T_len = 10
	p = 2
	q = 1

	# fix_receiver: tab = [theta, alpha_1, alpha_2], B = I
	true_alpha = c(0.5, 0.2)
	true_theta = 0.1
	true_tab   = c(true_theta, true_alpha)

	W = array(rnorm(m * m * p, sd = 0.3), dim = c(m, m, p))
	Z = array(rnorm(m * m * q * T_len, sd = 0.3), dim = c(m, m, q, T_len))
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
	ETA_t = sir::eta_tab(true_tab, W, X[,,t, drop=FALSE], Z[,,,t, drop=FALSE],
							fix_receiver = TRUE)
	lambda_t = exp(ETA_t[,,1])
	diag(lambda_t) = 0
	lambda_t = pmin(lambda_t, 50)
	lambda_t = pmax(lambda_t, 0.01)
	Y_t = matrix(rpois(m * m, lambda = c(lambda_t)), m, m)
	diag(Y_t) = NA
	Y[,,t] = Y_t
	}

	X_fit = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) {
	X_fit[,,t] = Y[,,t-1]
	X_fit[is.na(X_fit[,,t])] = 0
	}

	model = sir(Y = Y, W = W, X = X_fit, Z = Z,
				 family = "poisson", fix_receiver = TRUE,
				 calc_se = FALSE, trace = FALSE, max_iter = 50)

	est = unname(model$tab)

	# fix_receiver tab = [theta, alpha_1, alpha_2]
	expect_equal(est[1], true_theta, tolerance = 0.5,
				 label = "fix_receiver theta recovery")
	expect_equal(est[2], true_alpha[1], tolerance = 0.5,
				 label = "fix_receiver alpha_1 recovery")
	expect_equal(est[3], true_alpha[2], tolerance = 0.5,
				 label = "fix_receiver alpha_2 recovery")
})
