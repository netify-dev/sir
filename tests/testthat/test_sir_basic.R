test_that("SIR works with complete data", {
	set.seed(123)
	m = 10
	T_len = 5
	p = 2
	q = 2
	
	# generate synthetic data
	Y = array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
	W = array(rnorm(m * m * p), dim = c(m, m, p))
	X = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
	
	# test with Poisson family
	model_pois = sir(Y = Y, W = W, X = X, Z = Z, 
					family = "poisson", method = "ALS", 
					calc_se = FALSE, trace = FALSE, max_iter = 5)
	
	expect_s3_class(model_pois, "sir")
	expect_true(!is.null(model_pois$A))
	expect_true(!is.null(model_pois$B))
	expect_equal(dim(model_pois$A), c(m, m))
	expect_equal(dim(model_pois$B), c(m, m))
	
	# test with Normal family
	Y_normal = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	model_norm = sir(Y = Y_normal, W = W, X = X, Z = Z, 
					family = "normal", method = "ALS", 
					calc_se = FALSE, trace = FALSE, max_iter = 5)
	
	expect_s3_class(model_norm, "sir")
	expect_true(!is.null(model_norm$sigma2))
	
	# test with Binomial family
	Y_binom = array(rbinom(m * m * T_len, 1, 0.5), dim = c(m, m, T_len))
	model_binom = sir(Y = Y_binom, W = W, X = X, Z = Z, 
					 family = "binomial", method = "ALS", 
					 calc_se = FALSE, trace = FALSE, max_iter = 5)
	
	expect_s3_class(model_binom, "sir")
})

test_that("SIR works without influence matrices (W = NULL)", {
	set.seed(123)
	m = 10
	T_len = 5
	q = 2
	
	# generate data without W
	Y = array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
	
	# should work with W = NULL (no bilinear part)
	model = sir(Y = Y, W = NULL, X = NULL, Z = Z, 
				 family = "poisson", method = "ALS", 
				 calc_se = FALSE, trace = FALSE, max_iter = 5)
	
	expect_s3_class(model, "sir")
	expect_equal(dim(model$A), c(m, m))
	expect_equal(dim(model$B), c(m, m))
	# a and B should be zero matrices when W is NULL
	expect_true(all(model$A == 0))
	expect_true(all(model$B == 0))
	
	# check that theta is estimated
	expect_true(length(model$summ$coef) == q)
})

test_that("SIR works without exogenous covariates (Z = NULL)", {
	set.seed(123)
	m = 10
	T_len = 5
	p = 2
	
	# generate data without Z
	Y = array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
	W = array(rnorm(m * m * p), dim = c(m, m, p))
	X = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	
	# should work with Z = NULL
	model = sir(Y = Y, W = W, X = X, Z = NULL, 
				 family = "poisson", method = "ALS", 
				 calc_se = FALSE, trace = FALSE, max_iter = 5)
	
	expect_s3_class(model, "sir")
	expect_true(!is.null(model$A))
	expect_true(!is.null(model$B))
	
	# check that only alpha and beta are estimated (no theta)
	# p=2, so we should have p-1 alpha parameters + p beta parameters = 2p-1 = 3
	expect_true(length(model$summ$coef) == 2*p - 1)
})

test_that("SIR works with neither W nor Z (intercept-only model)", {
	set.seed(123)
	m = 10
	T_len = 5
	
	# generate data without W or Z
	Y = array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
	
	# should work with both W = NULL and Z = NULL
	model = sir(Y = Y, W = NULL, X = NULL, Z = NULL, 
				 family = "poisson", method = "ALS", 
				 calc_se = FALSE, trace = FALSE, max_iter = 5)
	
	expect_s3_class(model, "sir")
	# should have empty parameter vector
	expect_true(length(model$summ$coef) == 0)
})

test_that("SIR handles missing values (NA) in Y", {
	set.seed(123)
	m = 10
	T_len = 5
	p = 2
	q = 2
	
	# generate data with some NAs
	Y = array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
	# add some missing values
	Y[sample(1:length(Y), size = 20)] = NA
	
	W = array(rnorm(m * m * p), dim = c(m, m, p))
	X = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
	
	# should handle NAs gracefully
	model = sir(Y = Y, W = W, X = X, Z = Z, 
				 family = "poisson", method = "ALS", 
				 calc_se = FALSE, trace = FALSE, max_iter = 5)
	
	expect_s3_class(model, "sir")
	expect_true(!is.null(model$A))
	expect_true(!is.null(model$B))
})

test_that("SIR handles different fitting methods", {
	set.seed(123)
	m = 8
	T_len = 3
	p = 2
	q = 1
	
	# generate small dataset for faster testing
	Y = array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
	W = array(rnorm(m * m * p), dim = c(m, m, p))
	X = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
	
	# test ALS method
	model_als = sir(Y = Y, W = W, X = X, Z = Z, 
					 family = "poisson", method = "ALS", 
					 calc_se = FALSE, trace = FALSE, max_iter = 5)
	
	expect_s3_class(model_als, "sir")
	expect_equal(model_als$method, "ALS")
	
	# test optim method
	model_optim = sir(Y = Y, W = W, X = X, Z = Z, 
					 family = "poisson", method = "optim", 
					 calc_se = FALSE, trace = 0)
	
	expect_s3_class(model_optim, "sir")
	expect_equal(model_optim$method, "optim")
})

test_that("SIR computes standard errors when requested", {
	set.seed(123)
	m = 8
	T_len = 3
	p = 2
	q = 1
	
	# generate small dataset
	Y = array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
	W = array(rnorm(m * m * p), dim = c(m, m, p))
	X = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
	
	# test with SE calculation
	model_se = sir(Y = Y, W = W, X = X, Z = Z, 
					family = "poisson", method = "ALS", 
					calc_se = TRUE, trace = FALSE, max_iter = 5)
	
	expect_s3_class(model_se, "sir")
	expect_true("se" %in% colnames(model_se$summ))
	expect_true("rse" %in% colnames(model_se$summ))
	expect_true(all(!is.na(model_se$summ$coef)))
	
	# test without SE calculation
	model_no_se = sir(Y = Y, W = W, X = X, Z = Z, 
					 family = "poisson", method = "ALS", 
					 calc_se = FALSE, trace = FALSE, max_iter = 5)
	
	expect_s3_class(model_no_se, "sir")
	expect_false("se" %in% colnames(model_no_se$summ))
})

test_that("SIR handles 3D Z array correctly", {
	set.seed(123)
	m = 8
	T_len = 3
	p = 2
	
	Y = array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
	W = array(rnorm(m * m * p), dim = c(m, m, p))
	X = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	
	# test with 3D Z (should be treated as q=1)
	Z_3d = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	
	model = sir(Y = Y, W = W, X = X, Z = Z_3d, 
				 family = "poisson", method = "ALS", 
				 calc_se = FALSE, trace = FALSE, max_iter = 5)
	
	expect_s3_class(model, "sir")
	# should have 1 theta parameter
	expect_true(sum(grepl("^\\(Z\\)", rownames(model$summ))) == 1)
})

test_that("Print method works", {
	set.seed(123)
	m = 8
	T_len = 3
	p = 2
	q = 1
	
	Y = array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
	W = array(rnorm(m * m * p), dim = c(m, m, p))
	X = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
	
	model = sir(Y = Y, W = W, X = X, Z = Z, 
				 family = "poisson", method = "ALS", 
				 calc_se = FALSE, trace = FALSE, max_iter = 5)
	
	# should not error when printing (capture output to test it works)
	output = capture.output(print(model))
	expect_true(length(output) > 0)
})