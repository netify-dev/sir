test_that("rel_covar creates all three effects", {
	set.seed(42)
	m = 6
	T_len = 4
	arr = array(rnorm(m * m * T_len), dim = c(m, m, T_len))

	Z = rel_covar(arr, "trade")

	expect_equal(dim(Z), c(m, m, 3, T_len))
	expect_equal(dimnames(Z)[[3]], c("trade", "trade_recip", "trade_trans"))
})

test_that("rel_covar main effect matches input", {
	set.seed(42)
	m = 5
	T_len = 3
	arr = array(rnorm(m * m * T_len), dim = c(m, m, T_len))

	Z = rel_covar(arr, "x", effects = "main")

	expect_equal(dim(Z), c(m, m, 1, T_len))
	# main effect should equal the input
	for (t in 1:T_len) {
	expect_equal(Z[,,1,t], arr[,,t])
	}
})

test_that("rel_covar reciprocal is transpose", {
	set.seed(42)
	m = 5
	T_len = 3
	arr = array(rnorm(m * m * T_len), dim = c(m, m, T_len))

	Z = rel_covar(arr, "x", effects = "reciprocal")

	for (t in 1:T_len) {
	expect_equal(Z[,,1,t], t(arr[,,t]))
	}
})

test_that("rel_covar subset of effects works", {
	set.seed(42)
	arr = array(rnorm(5 * 5 * 3), dim = c(5, 5, 3))

	Z = rel_covar(arr, "y", effects = c("main", "transitive"))

	expect_equal(dim(Z)[3], 2)
	expect_equal(dimnames(Z)[[3]], c("y", "y_trans"))
})

test_that("rel_covar output works with sir()", {
	set.seed(42)
	m = 6
	T_len = 5
	Y = array(rpois(m * m * T_len, 2), dim = c(m, m, T_len))
	arr = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	W = array(rnorm(m * m * 2), dim = c(m, m, 2))
	X = array(0, dim = c(m, m, T_len))
	for (t in 2:T_len) X[,,t] = log(Y[,,t-1] + 1)

	Z = rel_covar(arr, "covar", effects = c("main", "reciprocal"))

	fit = sir(Y, W, X, Z = Z, family = "poisson",
			 fix_receiver = TRUE, calc_se = FALSE, max_iter = 5)

	expect_s3_class(fit, "sir")
	expect_equal(fit$q, 2)
})
