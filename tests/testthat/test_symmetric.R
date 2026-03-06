test_that("SIR handles symmetric (undirected) networks", {
	set.seed(42)
	m = 10
	T_len = 5
	p = 2
	q = 1

	# create symmetric Y
	Y = array(0, dim = c(m, m, T_len))
	for (t in 1:T_len) {
	Y_t = matrix(rpois(m * m, 2), m, m)
	Y[,,t] = (Y_t + t(Y_t)) / 2
	diag(Y[,,t]) = NA
	}

	W = array(rnorm(m * m * p), dim = c(m, m, p))
	X = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) X[,,t] = Y[,,t-1]
	X[is.na(X)] = 0

	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))

	model = sir(Y = Y, W = W, X = X, Z = Z,
				 family = "poisson", symmetric = TRUE,
				 calc_se = TRUE, trace = FALSE, max_iter = 10)

	expect_s3_class(model, "sir")
	expect_true(model$symmetric)
	expect_true(model$fix_receiver)

	# fitted values should be symmetric
	fv = fitted(model)
	for (t in 1:T_len) {
	fv_t = fv[,,t]
	fv_t[is.na(fv_t)] = 0
	expect_equal(fv_t, t(fv_t), tolerance = 1e-10)
	}

	# print and summary should work without error
	expect_no_error(capture.output(print(model)))
	s = summary(model)
	expect_true(s$symmetric)
})

test_that("symmetric warns when Y is not symmetric", {
	set.seed(42)
	m = 10
	T_len = 5

	# asymmetric Y
	Y = array(rpois(m * m * T_len, 2), dim = c(m, m, T_len))
	X = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) X[,,t] = Y[,,t-1]
	W = array(rnorm(m * m * 2), dim = c(m, m, 2))

	expect_message(
	sir(Y = Y, W = W, X = X, family = "poisson",
		symmetric = TRUE, calc_se = FALSE, max_iter = 5),
	"Symmetrizing Y"
	)
})

test_that("symmetric is incompatible with bipartite", {
	Y = array(1, dim = c(5, 6, 3))
	expect_error(
	sir(Y = Y, family = "poisson", symmetric = TRUE),
	"not compatible with bipartite"
	)
})
