test_that("SIR works with empty W array (p=0)", {
  set.seed(123)
  m <- 10
  T_len <- 5
  q <- 2
  
  # Generate data with empty W (p=0)
  Y <- array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
  W <- array(dim = c(m, m, 0))  # Empty W with p=0
  X <- array(rnorm(m * m * T_len), dim = c(m, m, T_len))
  Z <- array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
  
  # Should work with empty W
  model <- sir(Y = Y, W = W, X = X, Z = Z, 
               family = "poisson", method = "ALS", 
               calcSE = FALSE, trace = FALSE, max_iter = 5)
  
  expect_s3_class(model, "sir")
  expect_equal(dim(model$A), c(m, m))
  expect_equal(dim(model$B), c(m, m))
  # A and B should be zero matrices when p=0
  expect_true(all(model$A == 0))
  expect_true(all(model$B == 0))
  
  # Check that only theta is estimated
  expect_true(length(model$summ$coef) == q)
})

test_that("SIR handles empty W with optim method", {
  set.seed(123)
  m <- 8
  T_len <- 3
  q <- 1
  
  Y <- array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
  W <- array(numeric(0), dim = c(m, m, 0))  # Another way to create empty W
  X <- array(rnorm(m * m * T_len), dim = c(m, m, T_len))
  Z <- array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
  
  model <- sir(Y = Y, W = W, X = X, Z = Z, 
               family = "poisson", method = "optim", 
               calcSE = FALSE, trace = 0)
  
  expect_s3_class(model, "sir")
  expect_true(all(model$A == 0))
  expect_true(all(model$B == 0))
})

test_that("SIR handles empty W with all distribution families", {
  set.seed(123)
  m <- 8
  T_len <- 3
  q <- 2
  
  W <- array(dim = c(m, m, 0))
  X <- array(rnorm(m * m * T_len), dim = c(m, m, T_len))
  Z <- array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
  
  # Test Poisson
  Y_pois <- array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
  model_pois <- sir(Y = Y_pois, W = W, X = X, Z = Z, 
                    family = "poisson", method = "ALS", 
                    calcSE = FALSE, trace = FALSE, max_iter = 5)
  expect_s3_class(model_pois, "sir")
  
  # Test Normal
  Y_norm <- array(rnorm(m * m * T_len), dim = c(m, m, T_len))
  model_norm <- sir(Y = Y_norm, W = W, X = X, Z = Z, 
                    family = "normal", method = "ALS", 
                    calcSE = FALSE, trace = FALSE, max_iter = 5)
  expect_s3_class(model_norm, "sir")
  
  # Test Binomial
  Y_binom <- array(rbinom(m * m * T_len, 1, 0.5), dim = c(m, m, T_len))
  model_binom <- sir(Y = Y_binom, W = W, X = X, Z = Z, 
                     family = "binomial", method = "ALS", 
                     calcSE = FALSE, trace = FALSE, max_iter = 5)
  expect_s3_class(model_binom, "sir")
})

test_that("SIR handles empty W with missing data in Y", {
  set.seed(123)
  m <- 8
  T_len <- 3
  q <- 2
  
  Y <- array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
  Y[sample(1:length(Y), 10)] <- NA  # Add missing values
  
  W <- array(dim = c(m, m, 0))
  X <- array(rnorm(m * m * T_len), dim = c(m, m, T_len))
  Z <- array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
  
  model <- sir(Y = Y, W = W, X = X, Z = Z, 
               family = "poisson", method = "ALS", 
               calcSE = FALSE, trace = FALSE, max_iter = 5)
  
  expect_s3_class(model, "sir")
  expect_true(all(model$A == 0))
  expect_true(all(model$B == 0))
})

test_that("SIR handles empty W and NULL Z", {
  set.seed(123)
  m <- 8
  T_len <- 3
  
  Y <- array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
  W <- array(dim = c(m, m, 0))
  X <- array(rnorm(m * m * T_len), dim = c(m, m, T_len))
  
  # Both W empty and Z NULL - essentially an intercept-only model
  model <- sir(Y = Y, W = W, X = X, Z = NULL, 
               family = "poisson", method = "ALS", 
               calcSE = FALSE, trace = FALSE, max_iter = 5)
  
  expect_s3_class(model, "sir")
  expect_true(length(model$summ$coef) == 0)  # No parameters
})

test_that("SIR with empty W computes standard errors", {
  set.seed(123)
  m <- 6
  T_len <- 3
  q <- 1
  
  Y <- array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
  W <- array(dim = c(m, m, 0))
  X <- array(rnorm(m * m * T_len), dim = c(m, m, T_len))
  Z <- array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
  
  model <- sir(Y = Y, W = W, X = X, Z = Z, 
               family = "poisson", method = "ALS", 
               calcSE = TRUE, trace = FALSE, max_iter = 5)
  
  expect_s3_class(model, "sir")
  expect_true("se" %in% colnames(model$summ))
  expect_true("rse" %in% colnames(model$summ))
})