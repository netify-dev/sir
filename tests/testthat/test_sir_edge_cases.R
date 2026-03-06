test_that("SIR handles small networks", {
	set.seed(123)
	m = 3  # Very small network
	T_len = 2
	p = 1
	q = 1
	
	Y = array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
	W = array(rnorm(m * m * p), dim = c(m, m, p))
	X = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
	
	model = sir(Y = Y, W = W, X = X, Z = Z, 
				 family = "poisson", method = "ALS", 
				 calc_se = FALSE, trace = FALSE, max_iter = 5)
	
	expect_s3_class(model, "sir")
})

test_that("SIR handles single time point", {
	set.seed(123)
	m = 5
	T_len = 1  # Single time point
	p = 2
	q = 1

	Y = array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
	W = array(rnorm(m * m * p), dim = c(m, m, p))
	X = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))

	model = sir(Y = Y, W = W, X = X, Z = Z,
				 family = "poisson", method = "ALS",
				 calc_se = FALSE, trace = FALSE, max_iter = 5)

	expect_s3_class(model, "sir")
})

test_that("SIR warns with T=1 and zero X", {
	set.seed(123)
	m = 5
	T_len = 1
	p = 2
	q = 1

	Y = array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
	W = array(rnorm(m * m * p), dim = c(m, m, p))
	X = array(0, dim = c(m, m, T_len))
	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))

	expect_warning(
	sir(Y = Y, W = W, X = X, Z = Z,
		family = "poisson", method = "ALS",
		calc_se = FALSE, trace = FALSE, max_iter = 5),
	"Only one time period"
	)
})

test_that("SIR handles high proportion of missing data", {
	set.seed(123)
	m = 10
	T_len = 5
	p = 2
	q = 2
	
	Y = array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
	# add 50% missing values
	Y[sample(1:length(Y), size = length(Y)/2)] = NA
	
	W = array(rnorm(m * m * p), dim = c(m, m, p))
	X = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
	
	model = sir(Y = Y, W = W, X = X, Z = Z, 
				 family = "poisson", method = "ALS", 
				 calc_se = FALSE, trace = FALSE, max_iter = 5)
	
	expect_s3_class(model, "sir")
})

test_that("SIR handles diagonal missing pattern", {
	set.seed(123)
	m = 10
	T_len = 5
	p = 2
	q = 2
	
	Y = array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
	# set diagonal to NA (self-loops missing)
	for(t in 1:T_len) {
	diag(Y[,,t]) = NA
	}
	
	W = array(rnorm(m * m * p), dim = c(m, m, p))
	X = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
	
	model = sir(Y = Y, W = W, X = X, Z = Z, 
				 family = "poisson", method = "ALS", 
				 calc_se = FALSE, trace = FALSE, max_iter = 5)
	
	expect_s3_class(model, "sir")
	# check diagonal of A and B are zero (as enforced in the model)
	expect_true(all(diag(model$A) == 0))
	expect_true(all(diag(model$B) == 0))
})

test_that("SIR handles zero variance in covariates", {
	set.seed(123)
	m = 8
	T_len = 3
	p = 2
	q = 2
	
	Y = array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
	W = array(rnorm(m * m * p), dim = c(m, m, p))
	X = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	Z = array(0, dim = c(m, m, q, T_len))  # All zeros
	Z[,,1,] = 1  # First covariate is constant
	
	model = sir(Y = Y, W = W, X = X, Z = Z, 
				 family = "poisson", method = "ALS", 
				 calc_se = FALSE, trace = FALSE, max_iter = 5)
	
	expect_s3_class(model, "sir")
})

test_that("SIR handles sparse networks", {
	set.seed(123)
	m = 10
	T_len = 5
	p = 2
	q = 2
	
	# create sparse network (mostly zeros)
	Y = array(0, dim = c(m, m, T_len))
	# add a few edges
	for(t in 1:T_len) {
	Y[sample(1:m, 3), sample(1:m, 3), t] = rpois(9, lambda = 2)
	}
	
	W = array(rnorm(m * m * p), dim = c(m, m, p))
	X = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
	
	model = sir(Y = Y, W = W, X = X, Z = Z, 
				 family = "poisson", method = "ALS", 
				 calc_se = FALSE, trace = FALSE, max_iter = 5)
	
	expect_s3_class(model, "sir")
})

test_that("SIR handles extreme values in Y", {
	set.seed(123)
	m = 8
	T_len = 3
	p = 2
	q = 1
	
	# test with very large counts
	Y = array(rpois(m * m * T_len, lambda = 100), dim = c(m, m, T_len))
	W = array(rnorm(m * m * p), dim = c(m, m, p))
	X = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
	
	model = sir(Y = Y, W = W, X = X, Z = Z, 
				 family = "poisson", method = "ALS", 
				 calc_se = FALSE, trace = FALSE, max_iter = 5)
	
	expect_s3_class(model, "sir")
})

test_that("SIR rejects invalid inputs", {
	set.seed(123)
	m = 5
	T_len = 3
	
	Y = array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
	
	# test invalid family
	expect_error(
	sir(Y = Y, family = "invalid"),
	"family must be one of"
	)
	
	# test invalid method
	expect_error(
	sir(Y = Y, family = "poisson", method = "invalid"),
	"method must be"
	)
	
	# test non-square Y (bipartite - should work with fix_receiver)
	Y_bipartite = array(rpois(5 * 6 * 3, 2), dim = c(5, 6, 3))
	W_bp = array(rnorm(5 * 5 * 2), dim = c(5, 5, 2))
	X_bp = array(rnorm(5 * 6 * 3), dim = c(5, 6, 3))
	model_bp = sir(Y = Y_bipartite, W = W_bp, X = X_bp, family = "poisson",
					method = "ALS", calc_se = FALSE, max_iter = 5)
	expect_s3_class(model_bp, "sir")
	expect_true(model_bp$bipartite)
	
	# test dimension mismatch between Y and W
	W_invalid = array(1, dim = c(4, 4, 2))
	X = array(1, dim = c(m, m, T_len))
	expect_error(
	sir(Y = Y, W = W_invalid, X = X, family = "poisson"),
	"Dimension mismatch"
	)
	
	# test W without X
	W = array(1, dim = c(m, m, 2))
	expect_error(
	sir(Y = Y, W = W, X = NULL, family = "poisson"),
	"X must be provided"
	)
})

test_that("SIR handles single influence covariate (p=1)", {
	set.seed(123)
	m = 8
	T_len = 3
	p = 1  # Single influence covariate
	q = 2
	
	Y = array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
	W = array(rnorm(m * m * p), dim = c(m, m, p))
	X = array(rnorm(m * m * T_len), dim = c(m, m, T_len))
	Z = array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
	
	model = sir(Y = Y, W = W, X = X, Z = Z, 
				 family = "poisson", method = "ALS", 
				 calc_se = FALSE, trace = FALSE, max_iter = 5)
	
	expect_s3_class(model, "sir")
	# with p=1, we have alpha[1]=1 (fixed), so only beta[1] is estimated
	# plus q theta parameters = q + 1 total parameters
	expect_equal(length(model$summ$coef), q + p)
})