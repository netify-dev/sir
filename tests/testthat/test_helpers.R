test_that("eta_tab calculates linear predictor correctly", {
	set.seed(123)
	m = 5
	T_len = 3
	p = 2
	q = 2
	
	W = array(rnorm(m * m * p), dim = c(m, m, p))
	X = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
	
	# create parameter vector
	theta = rnorm(q)
	a = rnorm(p-1)  # alpha[-1]
	b = rnorm(p)     # beta
	tab = c(theta, a, b)
	
	eta = eta_tab(tab, W, X, Z)
	
	expect_equal(dim(eta), c(m, m, T_len))
	expect_false(any(is.na(eta)))
})

test_that("mll_sir calculates negative log-likelihood", {
	set.seed(123)
	m = 5
	T_len = 3
	p = 2
	q = 2
	
	Y = array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
	W = array(rnorm(m * m * p), dim = c(m, m, p))
	X = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
	
	# create parameter vector
	theta = rnorm(q)
	a = rnorm(p-1)
	b = rnorm(p)
	tab = c(theta, a, b)
	
	nll = mll_sir(tab, Y, W, X, Z, "poisson")
	
	expect_true(is.numeric(nll))
	expect_true(length(nll) == 1)
	expect_true(nll >= 0)  # NLL should be non-negative
	expect_false(is.na(nll))
})

test_that("flatten_Y works correctly", {
	m = 3
	T_len = 2
	
	Y = array(1:(m*m*T_len), dim = c(m, m, T_len))
	Y_flat = sir:::flatten_Y(Y)
	
	expect_equal(length(Y_flat), m * m * T_len)
	expect_equal(Y_flat, c(Y))
})

test_that("flatten_Z works correctly", {
	m = 3
	T_len = 2
	q = 2
	
	# test 4D array
	Z_4d = array(1:(m*m*q*T_len), dim = c(m, m, q, T_len))
	Z_flat = sir:::flatten_Z(Z_4d)
	
	expect_equal(dim(Z_flat), c(m*m*T_len, q))
	expect_true(ncol(Z_flat) == q)
	
	# test 3D array (treated as q=1)
	Z_3d = array(1:(m*m*T_len), dim = c(m, m, T_len))
	Z_flat_3d = sir:::flatten_Z(Z_3d)
	
	expect_equal(dim(Z_flat_3d), c(m*m*T_len, 1))
	
	# test NULL
	expect_null(sir:::flatten_Z(NULL))
})

test_that("prepare_Z_list works correctly", {
	m = 3
	T_len = 2
	q = 2
	
	# test 4D array
	Z_4d = array(rnorm(m*m*q*T_len), dim = c(m, m, q, T_len))
	Z_list = sir:::prepare_Z_list(Z_4d)
	
	expect_true(is.list(Z_list))
	expect_equal(length(Z_list), q)
	expect_equal(dim(Z_list[[1]]), c(m, m, T_len))
	
	# test 3D array
	Z_3d = array(rnorm(m*m*T_len), dim = c(m, m, T_len))
	Z_list_3d = sir:::prepare_Z_list(Z_3d)
	
	expect_equal(length(Z_list_3d), 1)
	expect_equal(dim(Z_list_3d[[1]]), c(m, m, T_len))
	
	# test NULL
	expect_equal(length(sir:::prepare_Z_list(NULL)), 0)
})

test_that("cast_array works correctly", {
	# create sample dyadic data
	dyad_data = expand.grid(i = 1:3, j = 1:3, t = 1:2)
	dyad_data$value = rpois(nrow(dyad_data), lambda = 2)
	
	arr = cast_array(dyad_data, "value")
	
	expect_equal(dim(arr), c(3, 3, 2))
	expect_false(any(is.na(arr)))
	
	# test monadic option
	dyad_data_monadic = dyad_data
	# set off-diagonal to 0
	dyad_data_monadic$value[dyad_data_monadic$i != dyad_data_monadic$j] = 0
	# set specific diagonal values for testing
	dyad_data_monadic$value[dyad_data_monadic$i == 1 & dyad_data_monadic$j == 1] = 5
	dyad_data_monadic$value[dyad_data_monadic$i == 2 & dyad_data_monadic$j == 2] = 3
	dyad_data_monadic$value[dyad_data_monadic$i == 3 & dyad_data_monadic$j == 3] = 4
	
	arr_monadic = cast_array(dyad_data_monadic, "value", monadic = TRUE, row = TRUE)
	# check that diagonal values are correctly set (unname to ignore names)
	expect_equal(unname(diag(arr_monadic[,,1])), c(5, 3, 4))
	expect_equal(unname(diag(arr_monadic[,,2])), c(5, 3, 4))
})

test_that("C++ helper functions work correctly", {
	m = 4
	T_len = 2
	p = 2
	
	# test cpp_amprod_W_v
	W = array(rnorm(m * m * p), dim = c(m, m, p))
	v = rnorm(p)
	
	result = sir:::cpp_amprod_W_v(W, v)
	expect_equal(dim(result), c(m, m))
	
	# test cpp_tprod_A_X_Bt
	X = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	A = matrix(rnorm(m * m), m, m)
	B = matrix(rnorm(m * m), m, m)
	
	result = sir:::cpp_tprod_A_X_Bt(X, A, B)
	expect_equal(dim(result), c(m, m, T_len))
	
	# test cpp_construct_Wbeta_design
	beta = rnorm(p)
	Wbeta_design = sir:::cpp_construct_Wbeta_design(W, X, beta)
	expect_equal(dim(Wbeta_design), c(m*m*T_len, p))
	
	# test cpp_construct_Walpha_design
	alpha = rnorm(p)
	Walpha_design = sir:::cpp_construct_Walpha_design(W, X, alpha)
	expect_equal(dim(Walpha_design), c(m*m*T_len, p))
})