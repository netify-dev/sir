test_that("SIR handles small networks", {
  set.seed(123)
  m <- 3  # Very small network
  T_len <- 2
  p <- 1
  q <- 1
  
  Y <- array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
  W <- array(rnorm(m * m * p), dim = c(m, m, p))
  X <- array(rnorm(m * m * T_len), dim = c(m, m, T_len))
  Z <- array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
  
  model <- sir(Y = Y, W = W, X = X, Z = Z, 
               family = "poisson", method = "ALS", 
               calcSE = FALSE, trace = FALSE, max_iter = 5)
  
  expect_s3_class(model, "sir")
})

test_that("SIR handles single time point", {
  set.seed(123)
  m <- 5
  T_len <- 1  # Single time point
  p <- 2
  q <- 1
  
  Y <- array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
  W <- array(rnorm(m * m * p), dim = c(m, m, p))
  X <- array(rnorm(m * m * T_len), dim = c(m, m, T_len))
  Z <- array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
  
  model <- sir(Y = Y, W = W, X = X, Z = Z, 
               family = "poisson", method = "ALS", 
               calcSE = FALSE, trace = FALSE, max_iter = 5)
  
  expect_s3_class(model, "sir")
})

test_that("SIR handles high proportion of missing data", {
  set.seed(123)
  m <- 10
  T_len <- 5
  p <- 2
  q <- 2
  
  Y <- array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
  # Add 50% missing values
  Y[sample(1:length(Y), size = length(Y)/2)] <- NA
  
  W <- array(rnorm(m * m * p), dim = c(m, m, p))
  X <- array(rnorm(m * m * T_len), dim = c(m, m, T_len))
  Z <- array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
  
  model <- sir(Y = Y, W = W, X = X, Z = Z, 
               family = "poisson", method = "ALS", 
               calcSE = FALSE, trace = FALSE, max_iter = 5)
  
  expect_s3_class(model, "sir")
})

test_that("SIR handles diagonal missing pattern", {
  set.seed(123)
  m <- 10
  T_len <- 5
  p <- 2
  q <- 2
  
  Y <- array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
  # Set diagonal to NA (self-loops missing)
  for(t in 1:T_len) {
    diag(Y[,,t]) <- NA
  }
  
  W <- array(rnorm(m * m * p), dim = c(m, m, p))
  X <- array(rnorm(m * m * T_len), dim = c(m, m, T_len))
  Z <- array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
  
  model <- sir(Y = Y, W = W, X = X, Z = Z, 
               family = "poisson", method = "ALS", 
               calcSE = FALSE, trace = FALSE, max_iter = 5)
  
  expect_s3_class(model, "sir")
  # Check diagonal of A and B are zero (as enforced in the model)
  expect_true(all(diag(model$A) == 0))
  expect_true(all(diag(model$B) == 0))
})

test_that("SIR handles zero variance in covariates", {
  set.seed(123)
  m <- 8
  T_len <- 3
  p <- 2
  q <- 2
  
  Y <- array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
  W <- array(rnorm(m * m * p), dim = c(m, m, p))
  X <- array(rnorm(m * m * T_len), dim = c(m, m, T_len))
  Z <- array(0, dim = c(m, m, q, T_len))  # All zeros
  Z[,,1,] <- 1  # First covariate is constant
  
  model <- sir(Y = Y, W = W, X = X, Z = Z, 
               family = "poisson", method = "ALS", 
               calcSE = FALSE, trace = FALSE, max_iter = 5)
  
  expect_s3_class(model, "sir")
})

test_that("SIR handles sparse networks", {
  set.seed(123)
  m <- 10
  T_len <- 5
  p <- 2
  q <- 2
  
  # Create sparse network (mostly zeros)
  Y <- array(0, dim = c(m, m, T_len))
  # Add a few edges
  for(t in 1:T_len) {
    Y[sample(1:m, 3), sample(1:m, 3), t] <- rpois(9, lambda = 2)
  }
  
  W <- array(rnorm(m * m * p), dim = c(m, m, p))
  X <- array(rnorm(m * m * T_len), dim = c(m, m, T_len))
  Z <- array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
  
  model <- sir(Y = Y, W = W, X = X, Z = Z, 
               family = "poisson", method = "ALS", 
               calcSE = FALSE, trace = FALSE, max_iter = 5)
  
  expect_s3_class(model, "sir")
})

test_that("SIR handles extreme values in Y", {
  set.seed(123)
  m <- 8
  T_len <- 3
  p <- 2
  q <- 1
  
  # Test with very large counts
  Y <- array(rpois(m * m * T_len, lambda = 100), dim = c(m, m, T_len))
  W <- array(rnorm(m * m * p), dim = c(m, m, p))
  X <- array(rnorm(m * m * T_len), dim = c(m, m, T_len))
  Z <- array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
  
  model <- sir(Y = Y, W = W, X = X, Z = Z, 
               family = "poisson", method = "ALS", 
               calcSE = FALSE, trace = FALSE, max_iter = 5)
  
  expect_s3_class(model, "sir")
})

test_that("SIR rejects invalid inputs", {
  set.seed(123)
  m <- 5
  T_len <- 3
  
  Y <- array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
  
  # Test invalid family
  expect_error(
    sir(Y = Y, family = "invalid"),
    "Invalid family"
  )
  
  # Test invalid method
  expect_error(
    sir(Y = Y, family = "poisson", method = "invalid"),
    "Invalid method"
  )
  
  # Test non-square Y
  Y_invalid <- array(1, dim = c(5, 6, 3))
  expect_error(
    sir(Y = Y_invalid, family = "poisson"),
    "Y must be a 3D array"
  )
  
  # Test dimension mismatch between Y and W
  W_invalid <- array(1, dim = c(4, 4, 2))
  X <- array(1, dim = c(m, m, T_len))
  expect_error(
    sir(Y = Y, W = W_invalid, X = X, family = "poisson"),
    "Dimensions.*must align"
  )
  
  # Test W without X
  W <- array(1, dim = c(m, m, 2))
  expect_error(
    sir(Y = Y, W = W, X = NULL, family = "poisson"),
    "X must be provided if W is provided"
  )
})

test_that("SIR handles single influence covariate (p=1)", {
  set.seed(123)
  m <- 8
  T_len <- 3
  p <- 1  # Single influence covariate
  q <- 2
  
  Y <- array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
  W <- array(rnorm(m * m * p), dim = c(m, m, p))
  X <- array(rnorm(m * m * T_len), dim = c(m, m, T_len))
  Z <- array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
  
  model <- sir(Y = Y, W = W, X = X, Z = Z, 
               family = "poisson", method = "ALS", 
               calcSE = FALSE, trace = FALSE, max_iter = 5)
  
  expect_s3_class(model, "sir")
  # With p=1, we have alpha[1]=1 (fixed), so only beta[1] is estimated
  # Plus q theta parameters = q + 1 total parameters
  expect_equal(length(model$summ$coef), q + p)
})