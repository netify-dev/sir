test_that("SIR handles symmetric (undirected) networks", {
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
	for (t in 2:T_len) X[,,t] = Y[,,t-1]
	X[is.na(X)] = 0

	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
	for (t in seq_len(T_len)) {
		for (k in seq_len(q)) {
			Z[, , k, t] = (Z[, , k, t] + t(Z[, , k, t])) / 2
		}
	}

	model = sir(Y = Y, W = W, X = X, Z = Z,
				 family = "poisson", symmetric = TRUE,
				 calc_se = TRUE, trace = FALSE, max_iter = 10)

	expect_s3_class(model, "sir")
	expect_true(model$symmetric)
	expect_true(model$fix_receiver)

		# check fitted-value symmetry
		fv = fitted(model)
		for (t in 1:T_len) {
		fv_t = fv[,,t]
		fv_t[is.na(fv_t)] = 0
		expect_equal(fv_t, t(fv_t), tolerance = 1e-10)
		}

		# newdata predictions should use the same mirrored symmetric convention
		pred_link = predict(model, newdata = list(W = W, X = X, Z = Z), type = "link")
		pred_resp = predict(model, newdata = list(W = W, X = X, Z = Z), type = "response")
		for (t in 1:T_len) {
			link_t = pred_link[, , t]
			resp_t = pred_resp[, , t]
			link_t[is.na(link_t)] = 0
			resp_t[is.na(resp_t)] = 0
			expect_equal(link_t, t(link_t), tolerance = 1e-10)
			expect_equal(resp_t, t(resp_t), tolerance = 1e-10)
			expect_true(all(is.na(diag(pred_link[, , t]))))
			expect_true(all(is.na(diag(pred_resp[, , t]))))
		}

		# check display methods
		expect_no_error(capture.output(print(model)))
	s = summary(model)
	expect_true(s$symmetric)
})

test_that("symmetric averages asymmetric normal outcomes", {
	set.seed(42)
	m = 10
	T_len = 5

	# build asymmetric continuous data
	Y = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	X = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) X[,,t] = Y[,,t-1]
	W = array(rnorm(m * m * 2), dim = c(m, m, 2))
	for (k in seq_len(2)) W[, , k] = (W[, , k] + t(W[, , k])) / 2

	expect_message(
	sir(Y = Y, W = W, X = X, family = "normal",
		symmetric = TRUE, calc_se = FALSE, max_iter = 5),
	"Symmetrizing Y"
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

	expect_error(
		sir(
			Y = Y,
			W = W,
			X = X,
			family = "poisson",
			symmetric = TRUE,
			calc_se = FALSE,
			max_iter = 5
		),
		"must already be symmetric"
	)
})

test_that("symmetric rejects asymmetric direct covariates", {
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
	W = array(0, dim = c(m, m, 1))
	W[, , 1] = 1
	diag(W[, , 1]) = 0
	Z = array(rnorm(m * m * T_len), dim = c(m, m, 1, T_len))

	expect_error(
		sir(
			Y = Y,
			W = W,
			X = X,
			Z = Z,
			family = "poisson",
			symmetric = TRUE,
			calc_se = FALSE,
			max_iter = 5
		),
		"direct-covariate"
	)
})

test_that("symmetric is incompatible with bipartite", {
	Y = array(1, dim = c(5, 6, 3))
	expect_error(
	sir(Y = Y, family = "poisson", symmetric = TRUE),
	"not compatible with bipartite"
	)
})
