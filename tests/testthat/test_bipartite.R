test_that("SIR handles bipartite networks with fix_receiver", {
	set.seed(42)
	n1 = 8   # senders
	n2 = 12  # receivers
	T_len = 5
	p = 2
	q = 1

	Y = array(rpois(n1 * n2 * T_len, 2), dim = c(n1, n2, T_len))
	W = array(rnorm(n1 * n1 * p), dim = c(n1, n1, p))  # sender covariates
	X = array(0, dim = c(n1, n2, T_len))
	for (t in 2:T_len) X[,,t] = log(Y[,,t-1] + 1)
	Z = array(rnorm(n1 * n2 * q * T_len), dim = c(n1, n2, q, T_len))

	model = sir(Y = Y, W = W, X = X, Z = Z,
				 family = "poisson", calc_se = TRUE, max_iter = 10)

	expect_s3_class(model, "sir")
	expect_true(model$bipartite)
	expect_true(model$fix_receiver)
	expect_equal(dim(model$A), c(n1, n1))
	expect_equal(dim(model$B), c(n2, n2))
	expect_equal(dim(fitted(model)), c(n1, n2, T_len))

	# print and summary should work without error
	expect_no_error(capture.output(print(model)))
	expect_no_error(summary(model))
})

test_that("Bipartite forces fix_receiver", {
	set.seed(42)
	n1 = 5
	n2 = 7
	T_len = 3

	Y = array(rpois(n1 * n2 * T_len, 2), dim = c(n1, n2, T_len))
	W = array(rnorm(n1 * n1 * 2), dim = c(n1, n1, 2))
	X = array(rnorm(n1 * n2 * T_len), dim = c(n1, n2, T_len))

	# should message about forcing fix_receiver
	expect_message(
	sir(Y = Y, W = W, X = X, family = "poisson",
		fix_receiver = FALSE, calc_se = FALSE, max_iter = 5),
	"fix_receiver"
	)
})

test_that("Bipartite rejects wrong W dimensions", {
	n1 = 5
	n2 = 7
	T_len = 3

	Y = array(1, dim = c(n1, n2, T_len))
	W_wrong = array(1, dim = c(n2, n2, 2))  # should be n1 x n1
	X = array(1, dim = c(n1, n2, T_len))

	expect_error(
	sir(Y = Y, W = W_wrong, X = X, family = "poisson"),
	"sender covariates"
	)
})

test_that("Bipartite works with normal family", {
	set.seed(42)
	n1 = 6
	n2 = 8
	T_len = 4

	Y = array(rnorm(n1 * n2 * T_len), dim = c(n1, n2, T_len))
	W = array(rnorm(n1 * n1 * 2), dim = c(n1, n1, 2))
	X = array(rnorm(n1 * n2 * T_len), dim = c(n1, n2, T_len))

	model = sir(Y = Y, W = W, X = X, family = "normal",
				 calc_se = FALSE, max_iter = 5)

	expect_s3_class(model, "sir")
	expect_true(model$bipartite)
})

test_that("Bipartite works without W", {
	set.seed(42)
	n1 = 5
	n2 = 7
	T_len = 3
	q = 1

	Y = array(rpois(n1 * n2 * T_len, 2), dim = c(n1, n2, T_len))
	Z = array(rnorm(n1 * n2 * q * T_len), dim = c(n1, n2, q, T_len))

	model = sir(Y = Y, W = NULL, X = NULL, Z = Z,
				 family = "poisson", calc_se = FALSE, max_iter = 5)

	expect_s3_class(model, "sir")
	expect_true(model$bipartite)
})
