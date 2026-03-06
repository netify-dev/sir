# test that C++ analytical gradient matches finite-difference gradient of mll_sir,
# and that diagonal handling is consistent between the R NLL and C++ gradient.

test_that("C++ gradient matches numerical gradient (Poisson)", {
	set.seed(42)
	m = 6
	T_len = 3
	p = 2
	q = 1

	Y = array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
	for (t in 1:T_len) diag(Y[,,t]) = NA
	W = array(rnorm(m * m * p, sd = 0.3), dim = c(m, m, p))
	X = array(rnorm(m * m * T_len, sd = 0.5), dim = c(m, m, T_len))
	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))

	# parameter vector: [theta, alpha_2:p, beta_1:p]
	tab = c(0.1, 0.2, 0.3, 0.1)

	Z_list = sir:::prepare_Z_list(Z)
	gH = sir:::cpp_mll_gH(tab, Y, W, X, Z_list, "poisson")
	analytical_grad = as.numeric(gH$grad)

	# numerical gradient via central differences
	eps = 1e-5
	numerical_grad = numeric(length(tab))
	for (i in seq_along(tab)) {
	tab_plus = tab
	tab_minus = tab
	tab_plus[i] = tab[i] + eps
	tab_minus[i] = tab[i] - eps
	nll_plus = mll_sir(tab_plus, Y, W, X, Z, "poisson")
	nll_minus = mll_sir(tab_minus, Y, W, X, Z, "poisson")
	numerical_grad[i] = (nll_plus - nll_minus) / (2 * eps)
	}

	expect_equal(analytical_grad, numerical_grad, tolerance = 1e-4,
				 label = "Poisson gradient consistency")
})

test_that("C++ gradient matches numerical gradient (Normal)", {
	set.seed(42)
	m = 6
	T_len = 3
	p = 2
	q = 1

	Y = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	for (t in 1:T_len) diag(Y[,,t]) = NA
	W = array(rnorm(m * m * p, sd = 0.3), dim = c(m, m, p))
	X = array(rnorm(m * m * T_len, sd = 0.5), dim = c(m, m, T_len))
	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))

	tab = c(0.1, 0.2, 0.3, 0.1)

	Z_list = sir:::prepare_Z_list(Z)
	gH = sir:::cpp_mll_gH(tab, Y, W, X, Z_list, "normal")
	analytical_grad = as.numeric(gH$grad)

	eps = 1e-5
	numerical_grad = numeric(length(tab))
	for (i in seq_along(tab)) {
	tab_plus = tab
	tab_minus = tab
	tab_plus[i] = tab[i] + eps
	tab_minus[i] = tab[i] - eps
	nll_plus = mll_sir(tab_plus, Y, W, X, Z, "normal")
	nll_minus = mll_sir(tab_minus, Y, W, X, Z, "normal")
	numerical_grad[i] = (nll_plus - nll_minus) / (2 * eps)
	}

	expect_equal(analytical_grad, numerical_grad, tolerance = 1e-4,
				 label = "Normal gradient consistency")
})

test_that("C++ gradient matches numerical gradient (Binomial)", {
	set.seed(42)
	m = 6
	T_len = 3
	p = 2
	q = 1

	Y = array(rbinom(m * m * T_len, 1, 0.4), dim = c(m, m, T_len))
	for (t in 1:T_len) diag(Y[,,t]) = NA
	W = array(rnorm(m * m * p, sd = 0.3), dim = c(m, m, p))
	X = array(rnorm(m * m * T_len, sd = 0.5), dim = c(m, m, T_len))
	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))

	tab = c(-0.1, 0.1, 0.2, -0.1)

	Z_list = sir:::prepare_Z_list(Z)
	gH = sir:::cpp_mll_gH(tab, Y, W, X, Z_list, "binomial")
	analytical_grad = as.numeric(gH$grad)

	eps = 1e-5
	numerical_grad = numeric(length(tab))
	for (i in seq_along(tab)) {
	tab_plus = tab
	tab_minus = tab
	tab_plus[i] = tab[i] + eps
	tab_minus[i] = tab[i] - eps
	nll_plus = mll_sir(tab_plus, Y, W, X, Z, "binomial")
	nll_minus = mll_sir(tab_minus, Y, W, X, Z, "binomial")
	numerical_grad[i] = (nll_plus - nll_minus) / (2 * eps)
	}

	expect_equal(analytical_grad, numerical_grad, tolerance = 1e-4,
				 label = "Binomial gradient consistency")
})

test_that("Diagonal exclusion is consistent between mll_sir and cpp_mll_gH", {
	set.seed(42)
	m = 5
	T_len = 2
	p = 2
	q = 1

	# y with NON-NA diagonals (the bug scenario)
	Y = array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
	# do NOT set diagonal to NA — test that both R and C++ handle this consistently
	W = array(rnorm(m * m * p, sd = 0.3), dim = c(m, m, p))
	X = array(rnorm(m * m * T_len, sd = 0.5), dim = c(m, m, T_len))
	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))

	tab = c(0.1, 0.2, 0.3, 0.1)

	Z_list = sir:::prepare_Z_list(Z)
	gH = sir:::cpp_mll_gH(tab, Y, W, X, Z_list, "poisson")
	analytical_grad = as.numeric(gH$grad)

	# mll_sir now also excludes diagonals, so finite differences should match
	eps = 1e-5
	numerical_grad = numeric(length(tab))
	for (i in seq_along(tab)) {
	tab_plus = tab
	tab_minus = tab
	tab_plus[i] = tab[i] + eps
	tab_minus[i] = tab[i] - eps
	nll_plus = mll_sir(tab_plus, Y, W, X, Z, "poisson")
	nll_minus = mll_sir(tab_minus, Y, W, X, Z, "poisson")
	numerical_grad[i] = (nll_plus - nll_minus) / (2 * eps)
	}

	expect_equal(analytical_grad, numerical_grad, tolerance = 1e-4,
				 label = "Diagonal-consistency gradient check")
})

test_that("C++ gradient matches numerical gradient with p=1", {
	set.seed(42)
	m = 6
	T_len = 3
	p = 1
	q = 2

	Y = array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
	for (t in 1:T_len) diag(Y[,,t]) = NA
	W = array(rnorm(m * m * p, sd = 0.3), dim = c(m, m, p))
	X = array(rnorm(m * m * T_len, sd = 0.5), dim = c(m, m, T_len))
	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))

	# p=1: tab = [theta_1, theta_2, beta_1] (alpha_1=1 fixed, no alpha_2:p)
	tab = c(0.1, -0.2, 0.3)

	Z_list = sir:::prepare_Z_list(Z)
	gH = sir:::cpp_mll_gH(tab, Y, W, X, Z_list, "poisson")
	analytical_grad = as.numeric(gH$grad)

	eps = 1e-5
	numerical_grad = numeric(length(tab))
	for (i in seq_along(tab)) {
	tab_plus = tab
	tab_minus = tab
	tab_plus[i] = tab[i] + eps
	tab_minus[i] = tab[i] - eps
	nll_plus = mll_sir(tab_plus, Y, W, X, Z, "poisson")
	nll_minus = mll_sir(tab_minus, Y, W, X, Z, "poisson")
	numerical_grad[i] = (nll_plus - nll_minus) / (2 * eps)
	}

	expect_equal(analytical_grad, numerical_grad, tolerance = 1e-4,
				 label = "p=1 gradient consistency")
})
