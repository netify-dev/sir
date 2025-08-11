test_that("eta_tab calculates linear predictor correctly", {
  set.seed(123)
  m <- 5
  T_len <- 3
  p <- 2
  q <- 2
  
  W <- array(rnorm(m * m * p), dim = c(m, m, p))
  X <- array(rnorm(m * m * T_len), dim = c(m, m, T_len))
  Z <- array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
  
  # Create parameter vector
  theta <- rnorm(q)
  a <- rnorm(p-1)  # alpha[-1]
  b <- rnorm(p)     # beta
  tab <- c(theta, a, b)
  
  eta <- eta_tab(tab, W, X, Z)
  
  expect_equal(dim(eta), c(m, m, T_len))
  expect_false(any(is.na(eta)))
})

test_that("mll_sir calculates negative log-likelihood", {
  set.seed(123)
  m <- 5
  T_len <- 3
  p <- 2
  q <- 2
  
  Y <- array(rpois(m * m * T_len, lambda = 2), dim = c(m, m, T_len))
  W <- array(rnorm(m * m * p), dim = c(m, m, p))
  X <- array(rnorm(m * m * T_len), dim = c(m, m, T_len))
  Z <- array(rnorm(m * m * q * T_len), dim = c(m, m, q, T_len))
  
  # Create parameter vector
  theta <- rnorm(q)
  a <- rnorm(p-1)
  b <- rnorm(p)
  tab <- c(theta, a, b)
  
  nll <- mll_sir(tab, Y, W, X, Z, "poisson")
  
  expect_true(is.numeric(nll))
  expect_true(length(nll) == 1)
  expect_true(nll >= 0)  # NLL should be non-negative
  expect_false(is.na(nll))
})

test_that("flatten_Y works correctly", {
  m <- 3
  T_len <- 2
  
  Y <- array(1:(m*m*T_len), dim = c(m, m, T_len))
  Y_flat <- flatten_Y(Y)
  
  expect_equal(length(Y_flat), m * m * T_len)
  expect_equal(Y_flat, c(Y))
})

test_that("flatten_Z works correctly", {
  m <- 3
  T_len <- 2
  q <- 2
  
  # Test 4D array
  Z_4d <- array(1:(m*m*q*T_len), dim = c(m, m, q, T_len))
  Z_flat <- flatten_Z(Z_4d)
  
  expect_equal(dim(Z_flat), c(m*m*T_len, q))
  expect_true(ncol(Z_flat) == q)
  
  # Test 3D array (treated as q=1)
  Z_3d <- array(1:(m*m*T_len), dim = c(m, m, T_len))
  Z_flat_3d <- flatten_Z(Z_3d)
  
  expect_equal(dim(Z_flat_3d), c(m*m*T_len, 1))
  
  # Test NULL
  expect_null(flatten_Z(NULL))
})

test_that("prepare_Z_list works correctly", {
  m <- 3
  T_len <- 2
  q <- 2
  
  # Test 4D array
  Z_4d <- array(rnorm(m*m*q*T_len), dim = c(m, m, q, T_len))
  Z_list <- prepare_Z_list(Z_4d)
  
  expect_true(is.list(Z_list))
  expect_equal(length(Z_list), q)
  expect_equal(dim(Z_list[[1]]), c(m, m, T_len))
  
  # Test 3D array
  Z_3d <- array(rnorm(m*m*T_len), dim = c(m, m, T_len))
  Z_list_3d <- prepare_Z_list(Z_3d)
  
  expect_equal(length(Z_list_3d), 1)
  expect_equal(dim(Z_list_3d[[1]]), c(m, m, T_len))
  
  # Test NULL
  expect_equal(length(prepare_Z_list(NULL)), 0)
})

test_that("castArray works correctly", {
  # Create sample dyadic data
  dyadData <- expand.grid(i = 1:3, j = 1:3, t = 1:2)
  dyadData$value <- rpois(nrow(dyadData), lambda = 2)
  
  arr <- castArray(dyadData, "value")
  
  expect_equal(dim(arr), c(3, 3, 2))
  expect_false(any(is.na(arr)))
  
  # Test monadic option
  dyadData_monadic <- dyadData
  # Set off-diagonal to 0
  dyadData_monadic$value[dyadData_monadic$i != dyadData_monadic$j] <- 0
  # Set specific diagonal values for testing
  dyadData_monadic$value[dyadData_monadic$i == 1 & dyadData_monadic$j == 1] <- 5
  dyadData_monadic$value[dyadData_monadic$i == 2 & dyadData_monadic$j == 2] <- 3
  dyadData_monadic$value[dyadData_monadic$i == 3 & dyadData_monadic$j == 3] <- 4
  
  arr_monadic <- castArray(dyadData_monadic, "value", monadic = TRUE, row = TRUE)
  # Check that diagonal values are correctly set (unname to ignore names)
  expect_equal(unname(diag(arr_monadic[,,1])), c(5, 3, 4))
  expect_equal(unname(diag(arr_monadic[,,2])), c(5, 3, 4))
})

test_that("C++ helper functions work correctly", {
  m <- 4
  T_len <- 2
  p <- 2
  
  # Test cpp_amprod_W_v
  W <- array(rnorm(m * m * p), dim = c(m, m, p))
  v <- rnorm(p)
  
  result <- cpp_amprod_W_v(W, v)
  expect_equal(dim(result), c(m, m))
  
  # Test cpp_tprod_A_X_Bt
  X <- array(rnorm(m * m * T_len), dim = c(m, m, T_len))
  A <- matrix(rnorm(m * m), m, m)
  B <- matrix(rnorm(m * m), m, m)
  
  result <- cpp_tprod_A_X_Bt(X, A, B)
  expect_equal(dim(result), c(m, m, T_len))
  
  # Test cpp_construct_Wbeta_design
  beta <- rnorm(p)
  Wbeta_design <- cpp_construct_Wbeta_design(W, X, beta)
  expect_equal(dim(Wbeta_design), c(m*m*T_len, p))
  
  # Test cpp_construct_Walpha_design
  alpha <- rnorm(p)
  Walpha_design <- cpp_construct_Walpha_design(W, X, alpha)
  expect_equal(dim(Walpha_design), c(m*m*T_len, p))
})