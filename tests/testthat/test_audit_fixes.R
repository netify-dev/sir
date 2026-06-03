test_that("als ignores one-mode square diagonals", {
	set.seed(101)
	n = 5
	t_len = 4
	Y = array(rnorm(n * n * t_len), dim = c(n, n, t_len))
	W = array(rnorm(n * n), dim = c(n, n, 1))
	X = array(rnorm(n * n * t_len), dim = c(n, n, t_len))

	fit_base = sir(
		Y = Y,
		W = W,
		X = X,
		family = "normal",
		calc_se = FALSE,
		seed = 1,
		max_iter = 20
	)

	Y_diag = Y
	W_diag = W
	X_diag = X
	for (t in seq_len(t_len)) {
		diag(Y_diag[, , t]) = 1000
		diag(X_diag[, , t]) = 1000
	}
	diag(W_diag[, , 1]) = 1000

	fit_diag = sir(
		Y = Y_diag,
		W = W_diag,
		X = X_diag,
		family = "normal",
		calc_se = FALSE,
		seed = 1,
		max_iter = 20
	)

	expect_equal(fit_diag$tab, fit_base$tab, tolerance = 1e-8)
	expect_true(all(is.na(diag(fit_diag$Y[, , 1]))))
	expect_true(all(diag(fit_diag$W[, , 1]) == 0))
	expect_true(all(diag(fit_diag$X[, , 1]) == 0))
})

test_that("rectangular 3d z is accepted for sender-side bipartite fits", {
	set.seed(102)
	n1 = 3
	n2 = 4
	t_len = 5
	Y = array(rnorm(n1 * n2 * t_len), dim = c(n1, n2, t_len))
	W = array(rnorm(n1 * n1), dim = c(n1, n1, 1))
	X = array(rnorm(n1 * n2 * t_len), dim = c(n1, n2, t_len))
	Z = array(rnorm(n1 * n2 * t_len), dim = c(n1, n2, t_len))

	fit = sir(
		Y = Y,
		W = W,
		X = X,
		Z = Z,
		family = "normal",
		calc_se = FALSE,
		seed = 1
	)

	expect_true(fit$bipartite)
	expect_equal(dim(fit$Z), c(n1, n2, 1, t_len))
	expect_equal(fit$nobs, n1 * n2 * t_len)
})

test_that("full-bilinear bipartite keeps square two-mode diagonals", {
	set.seed(103)
	n = 5
	t_len = 6
	Y = array(rnorm(n * n * t_len), dim = c(n, n, t_len))
	W = array(rnorm(n * n), dim = c(n, n, 1))
	W_recv = array(rnorm(n * n), dim = c(n, n, 1))
	X = array(rnorm(n * n * t_len), dim = c(n, n, t_len))

	fit = sir(
		Y = Y,
		W = W,
		W_recv = W_recv,
		X = X,
		family = "normal",
		calc_se = FALSE,
		seed = 1,
		max_iter = 20,
		n_restarts = 1
	)

	expect_equal(fit$nobs, n * n * t_len)
	expect_equal(sum(is.na(fit$Y)), 0)
})

test_that("full-bilinear bipartite validates and propagates Z missingness", {
	set.seed(105)
	n1 = 4
	n2 = 5
	t_len = 4
	Y = array(rnorm(n1 * n2 * t_len), dim = c(n1, n2, t_len))
	W = array(rnorm(n1 * n1), dim = c(n1, n1, 1))
	Wr = array(rnorm(n2 * n2), dim = c(n2, n2, 1))
	X = array(rnorm(n1 * n2 * t_len), dim = c(n1, n2, t_len))
	Z = array(rnorm(n1 * n2 * t_len), dim = c(n1, n2, t_len))
	Z[1, 2, 3] = NA

	fit = sir(
		Y = Y,
		W = W,
		W_recv = Wr,
		X = X,
		Z = Z,
		family = "normal",
		calc_se = FALSE,
		seed = 1,
		max_iter = 20,
		n_restarts = 1
	)

	expect_equal(dim(fit$Z), c(n1, n2, 1, t_len))
	expect_true(is.na(fit$Y[1, 2, 3]))
	expect_equal(fit$nobs, n1 * n2 * t_len - 1)
	expect_error(
		sir(Y = Y, W = W, W_recv = Wr, X = X[, , 1:3], family = "normal"),
		"X"
	)
})

test_that("full-bilinear bipartite aligns named sender and receiver arrays", {
	set.seed(106)
	send = paste0("s", 1:4)
	recv = paste0("r", 1:5)
	t_len = 5
	Y = array(rnorm(4 * 5 * t_len), dim = c(4, 5, t_len),
			  dimnames = list(send, recv, NULL))
	W0 = array(rnorm(4 * 4), dim = c(4, 4, 1),
			   dimnames = list(rev(send), rev(send), "w"))
	Wr0 = array(rnorm(5 * 5), dim = c(5, 5, 1),
				dimnames = list(rev(recv), rev(recv), "wr"))
	X0 = array(rnorm(4 * 5 * t_len), dim = c(4, 5, t_len),
			   dimnames = list(rev(send), rev(recv), NULL))

	fit = sir(Y, W = W0, W_recv = Wr0, X = X0, family = "normal",
			  calc_se = FALSE, seed = 1, max_iter = 20, n_restarts = 1)
	expect_equal(rownames(fit$A), send)
	expect_equal(colnames(fit$A), send)
	expect_equal(rownames(fit$B), recv)
	expect_equal(colnames(fit$B), recv)
})

test_that("normal optim loglik uses estimated gaussian scale", {
	set.seed(104)
	n = 5
	t_len = 5
	Z = array(rnorm(n * n * t_len), dim = c(n, n, 1, t_len))
	Y = array(0.6 * Z[, , 1, ] + rnorm(n * n * t_len, sd = 0.5), dim = c(n, n, t_len))

	fit = sir(
		Y = Y,
		Z = Z,
		family = "normal",
		method = "optim",
		calc_se = FALSE,
		seed = 1
	)

	rss = sum((fit$Y - fit$fitted.values)^2, na.rm = TRUE)
	expected_sigma2 = rss / (fit$nobs - length(fit$tab))
	expected_sigma2_mle = rss / fit$nobs
	expected_ll = sum(
		stats::dnorm(
			fit$Y,
			mean = fit$fitted.values,
			sd = sqrt(expected_sigma2_mle),
			log = TRUE
		),
		na.rm = TRUE
	)

	expect_equal(fit$sigma2, expected_sigma2, tolerance = 1e-10)
	expect_equal(fit$sigma2_mle, expected_sigma2_mle, tolerance = 1e-10)
	expect_equal(as.numeric(logLik(fit)), expected_ll, tolerance = 1e-10)
})
