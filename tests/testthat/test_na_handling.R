# tests for NA handling across Y, X, W, and Z inputs

test_that("sir handles NAs in Y (dyad-level missingness)", {
	set.seed(42)
	m = 10
	T_len = 5
	p = 2

	Y = array(rpois(m * m * T_len, 2), dim = c(m, m, T_len))
	for (t in 1:T_len) diag(Y[,,t]) = NA
	W = array(rnorm(m * m * p, sd = 0.2), dim = c(m, m, p))
	X = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) X[,,t] = Y[,,t-1]
	X[is.na(X)] = 0

	# randomly set 10% of off-diagonal entries to NA
	for (t in 1:T_len) {
	mask = matrix(runif(m * m) < 0.1, m, m)
	diag(mask) = FALSE
	Y[,,t][mask] = NA
	}

	expect_message({
	model = sir(Y = Y, W = W, X = X, family = "poisson",
				 fix_receiver = TRUE, calc_se = FALSE, trace = FALSE, max_iter = 10)
	}, "off-diagonal dyad-observations")
	expect_s3_class(model, "sir")
	expect_true(model$nobs < m * (m - 1) * T_len)
})

test_that("sir handles NAs in X (zeroed out, no corruption)", {
	set.seed(42)
	m = 10
	T_len = 5
	p = 2

	Y = array(rpois(m * m * T_len, 2), dim = c(m, m, T_len))
	for (t in 1:T_len) diag(Y[,,t]) = NA
	W = array(rnorm(m * m * p, sd = 0.2), dim = c(m, m, p))

	# x has off-diagonal NAs (e.g., from lagged Y with missing dyads)
	X = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) X[,,t] = Y[,,t-1]
	# intentionally leave NAs in X, and add extra off-diagonal NAs
	X[2, 5, 3] = NA
	X[7, 3, 4] = NA

	expect_message({
	model = sir(Y = Y, W = W, X = X, family = "poisson",
				 fix_receiver = TRUE, calc_se = FALSE, trace = FALSE, max_iter = 10)
	}, "Replacing.*NA values in X")
	expect_s3_class(model, "sir")
	# coefficients should be finite (NAs in X didn't corrupt)
	expect_true(all(is.finite(model$tab)))
})

test_that("sir handles NAs in W (zeroed out)", {
	set.seed(42)
	m = 10
	T_len = 5
	p = 2

	Y = array(rpois(m * m * T_len, 2), dim = c(m, m, T_len))
	for (t in 1:T_len) diag(Y[,,t]) = NA
	W = array(rnorm(m * m * p, sd = 0.2), dim = c(m, m, p))
	X = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) { X[,,t] = Y[,,t-1]; X[is.na(X[,,t])] = 0 }

	# set some W entries to NA (e.g., missing covariate for some dyads)
	W[1, 3, 1] = NA
	W[5, 7, 2] = NA

	expect_message({
	model = sir(Y = Y, W = W, X = X, family = "poisson",
				 fix_receiver = TRUE, calc_se = FALSE, trace = FALSE, max_iter = 10)
	}, "Replacing.*NA values in W")
	expect_s3_class(model, "sir")
	expect_true(all(is.finite(model$tab)))
})

test_that("sir handles NAs in Z (propagated to Y, then zeroed)", {
	set.seed(42)
	m = 10
	T_len = 5
	p = 2
	q = 1

	Y = array(rpois(m * m * T_len, 2), dim = c(m, m, T_len))
	for (t in 1:T_len) diag(Y[,,t]) = NA
	W = array(rnorm(m * m * p, sd = 0.2), dim = c(m, m, p))
	X = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) { X[,,t] = Y[,,t-1]; X[is.na(X[,,t])] = 0 }
	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))

	# set some Z entries to NA
	Z[2, 4, 1, 3] = NA
	Z[6, 8, 1, 5] = NA

	expect_message({
	model = sir(Y = Y, W = W, X = X, Z = Z, family = "poisson",
				 fix_receiver = TRUE, calc_se = FALSE, trace = FALSE, max_iter = 10)
	}, "Z contains NAs.*excluding.*dyad-observations")
	expect_s3_class(model, "sir")
	expect_true(all(is.finite(model$tab)))
})

test_that("sir handles listwise deletion (entire actor rows/cols NA)", {
	set.seed(42)
	m = 10
	T_len = 5
	p = 2

	Y = array(rpois(m * m * T_len, 2), dim = c(m, m, T_len))
	for (t in 1:T_len) diag(Y[,,t]) = NA
	W = array(rnorm(m * m * p, sd = 0.2), dim = c(m, m, p))
	X = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) { X[,,t] = Y[,,t-1]; X[is.na(X[,,t])] = 0 }

	# actor 3 is missing at time 4 (listwise deletion)
	Y[3, , 4] = NA
	Y[, 3, 4] = NA

	expect_message({
	model = sir(Y = Y, W = W, X = X, family = "poisson",
				 fix_receiver = TRUE, calc_se = FALSE, trace = FALSE, max_iter = 10)
	}, "off-diagonal dyad-observations")
	expect_s3_class(model, "sir")
	# fitted values at the NA positions should still be computed (eta is defined)
	# but Y was NA so those positions were excluded from the likelihood
	expect_true(all(is.finite(model$tab)))
})

test_that("sir with no NAs produces no NA messages", {
	set.seed(42)
	m = 8
	T_len = 3
	p = 2

	Y = array(rpois(m * m * T_len, 2), dim = c(m, m, T_len))
	for (t in 1:T_len) diag(Y[,,t]) = NA
	W = array(rnorm(m * m * p, sd = 0.2), dim = c(m, m, p))
	X = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) { X[,,t] = Y[,,t-1]; X[is.na(X[,,t])] = 0 }

	# no extra NAs beyond diagonal -- should produce no "missing" messages
	output = capture.output(
	model = sir(Y = Y, W = W, X = X, family = "poisson",
				 fix_receiver = TRUE, calc_se = FALSE, trace = FALSE, max_iter = 10),
	type = "message"
	)
	# should not contain "missing" or "NA" messages
	expect_false(any(grepl("missing|NA values|excluding", output)))
})

test_that("sir handles NAs in bipartite networks", {
	set.seed(42)
	n1 = 8
	n2 = 12
	T_len = 5
	p = 2

	Y = array(rpois(n1 * n2 * T_len, 2), dim = c(n1, n2, T_len))
	# set some entries to NA
	Y[1, 3, 2] = NA
	Y[5, , 4] = NA  # entire sender row missing at time 4
	W = array(rnorm(n1 * n1 * p, sd = 0.2), dim = c(n1, n1, p))
	X = array(0, dim = c(n1, n2, T_len))
	for (t in 2:T_len) X[,,t] = Y[,,t-1]
	X[is.na(X)] = 0

	expect_message({
	model = sir(Y = Y, W = W, X = X, family = "poisson",
				 fix_receiver = TRUE, calc_se = FALSE, trace = FALSE, max_iter = 10)
	}, "dyad-observations")
	expect_s3_class(model, "sir")
})
