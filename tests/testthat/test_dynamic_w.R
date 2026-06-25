test_that("dynamic W with a single covariate (p = 1) fits and gives valid SEs", {
	# regression: W[,,,t] collapses to a 2D matrix when p = 1, which used to crash
	# cpp_amprod_W_v ("Input array must have exactly 3 dimensions") inside the A/B
	# construction and eta_tab. it must stay 3D.
	set.seed(33)
	m = 10; T_len = 40
	W = array(0, dim = c(m, m, 1, T_len))
	W[, , 1, 1] = rnorm(m * m)
	for (t in 2:T_len) {
		W[, , 1, t] = 0.9 * W[, , 1, t - 1] + sqrt(0.19) * matrix(rnorm(m * m), m, m)
	}
	for (t in seq_len(T_len)) diag(W[, , 1, t]) = 0
	Y = array(0, c(m, m, T_len)); X = array(0, c(m, m, T_len))
	for (t in seq_len(T_len)) {
		if (t > 1) { lg = Y[, , t - 1]; lg[is.na(lg)] = 0; X[, , t] = lg / (m - 1) }
		Y[, , t] = 1.3 * W[, , 1, t] %*% X[, , t] + matrix(rnorm(m * m, sd = 0.4), m, m)
		diag(Y[, , t]) = NA
	}

	# fix_receiver p = 1: one free alpha; must not crash and must recover
	fit = suppressWarnings(sir(Y, W = W, X = X, family = "normal",
							   fix_receiver = TRUE, seed = 1))
	expect_equal(dim(fit$A), c(m, m, T_len))
	expect_equal(length(coef(fit)), 1L)
	expect_true(abs(coef(fit)[1] - 1.3) < 0.2)

	# eta_tab reconstructs the fitted values exactly (the A/B slice stayed 3D)
	et = sir:::eta_tab(fit$tab, fit$W, fit$X, fit$Z, fix_receiver = TRUE)
	off = vapply(seq_len(T_len), function(t) {
		d = abs(et[, , t] - fit$fitted.values[, , t]); diag(d) = NA
		max(d, na.rm = TRUE)
	}, numeric(1))
	expect_lt(max(off), 1e-8)

	# the SE machinery works: cluster sandwich is finite/PSD with a t(G-1) reference
	V = vcov(fit, type = "cluster")
	expect_true(all(is.finite(V)))
	expect_equal(attr(V, "cluster_df"), m - 1L)
	expect_true(all(is.finite(suppressMessages(confint(fit)))))

	# directed p = 1 (alpha fixed at 1, one free beta) also fits
	Z = array(rnorm(m * m * T_len), c(m, m, T_len))
	set.seed(7); Xd = array(rnorm(m * m * T_len), c(m, m, T_len)) / (m - 1)
	Yd = array(0, c(m, m, T_len))
	for (t in seq_len(T_len)) {
		Yd[, , t] = 0.5 * Z[, , t] +
			W[, , 1, t] %*% Xd[, , t] %*% t(1.4 * W[, , 1, t]) +
			matrix(rnorm(m * m, sd = 0.3), m, m)
		diag(Yd[, , t]) = NA
	}
	fd = suppressWarnings(sir(Yd, W = W, X = Xd, Z = Z, family = "normal", seed = 1))
	expect_equal(length(coef(fd)), 2L)
	expect_true(all(is.finite(vcov(fd, type = "cluster"))))
})

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
