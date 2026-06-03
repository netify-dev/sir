test_that("sir() accepts 4D W (dynamic influence covariates)", {
	set.seed(42)
	m = 6
	T_len = 5
	p = 2
	q = 1

	Y = array(rpois(m * m * T_len, 2), dim = c(m, m, T_len))
	# 4D W: influence covariates change over time
	W = array(rnorm(m * m * p * T_len), dim = c(m, m, p, T_len))
	X = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) X[,,t] = log(Y[,,t-1] + 1)
	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))

	fit = sir(Y, W, X, Z, family = "poisson",
			 fix_receiver = TRUE, calc_se = FALSE, max_iter = 5)

	expect_s3_class(fit, "sir")
	expect_true(fit$dynamic_W)
	# a should be 3D (m x m x T) for dynamic W
	expect_equal(length(dim(fit$A)), 3)
	expect_equal(dim(fit$A), c(m, m, T_len))
})

test_that("dynamic W with full bilinear (fix_receiver = FALSE)", {
	set.seed(42)
	m = 6
	T_len = 5
	p = 2

	Y = array(rpois(m * m * T_len, 2), dim = c(m, m, T_len))
	W = array(rnorm(m * m * p * T_len), dim = c(m, m, p, T_len))
	X = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) X[,,t] = log(Y[,,t-1] + 1)

	expect_warning({
		fit = sir(Y, W, X, family = "poisson",
				 fix_receiver = FALSE, calc_se = FALSE, max_iter = 5)
	}, "ALS did not converge")

	expect_s3_class(fit, "sir")
	expect_true(fit$dynamic_W)
	# both A and B should be 3D
	expect_equal(length(dim(fit$A)), 3)
	expect_equal(length(dim(fit$B)), 3)
})

test_that("dynamic W prediction works", {
	set.seed(42)
	m = 6
	T_len = 5
	p = 2

	Y = array(rpois(m * m * T_len, 2), dim = c(m, m, T_len))
	W = array(rnorm(m * m * p * T_len), dim = c(m, m, p, T_len))
	X = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) X[,,t] = log(Y[,,t-1] + 1)

	fit = sir(Y, W, X, family = "poisson",
			 fix_receiver = TRUE, calc_se = FALSE, max_iter = 5)

	pred = predict(fit)
	expect_equal(dim(pred), c(m, m, T_len))
	expect_true(all(pred >= 0, na.rm = TRUE))
})

test_that("rectangular bipartite dynamic W uses sender-by-sender dimensions", {
	set.seed(43)
	n1 = 5
	n2 = 7
	T_len = 4
	p = 2

	Y = array(rpois(n1 * n2 * T_len, 2), dim = c(n1, n2, T_len))
	W = array(rnorm(n1 * n1 * p * T_len), dim = c(n1, n1, p, T_len))
	X = array(rnorm(n1 * n2 * T_len), dim = c(n1, n2, T_len))

	fit = sir(Y, W, X, family = "poisson",
			  fix_receiver = FALSE, calc_se = FALSE, max_iter = 5)
	expect_true(fit$bipartite)
	expect_true(fit$dynamic_W)
	expect_equal(dim(fit$A), c(n1, n1, T_len))
	expect_equal(dim(fit$B), c(n2, n2, T_len))
})

test_that("dynamic W forces ALS method", {
	set.seed(42)
	m = 5
	T_len = 3
	p = 2

	Y = array(rpois(m * m * T_len, 2), dim = c(m, m, T_len))
	W = array(rnorm(m * m * p * T_len), dim = c(m, m, p, T_len))
	X = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) X[,,t] = log(Y[,,t-1] + 1)

	# dynamic W with method="optim" should switch to ALS
	suppressWarnings(
		expect_message({
		fit = sir(Y, W, X, family = "poisson",
					 method = "optim", fix_receiver = FALSE,
					 calc_se = FALSE, max_iter = 3)
		}, "ALS")
	)
	expect_equal(fit$method, "ALS")
})

test_that("summary and plot work with dynamic W", {
	set.seed(42)
	m = 6
	T_len = 5
	p = 2

	Y = array(rpois(m * m * T_len, 2), dim = c(m, m, T_len))
	W = array(rnorm(m * m * p * T_len), dim = c(m, m, p, T_len))
	X = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) X[,,t] = log(Y[,,t-1] + 1)

	fit = sir(Y, W, X, family = "poisson",
			 fix_receiver = TRUE, calc_se = FALSE, max_iter = 5)

	# summary should not error
	s = summary(fit)
	expect_s3_class(s, "summary.sir")

	# print should not error
	output = capture.output(print(fit))
	expect_true(length(output) > 0)

	expect_no_error(plot(fit, which = 1:4, combine = FALSE))
})
