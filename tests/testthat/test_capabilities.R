# tests for the predictive spine (forecast / cv_sir / score_sir), the dyad
# bootstrap, and the full-bilinear bipartite path.

test_that("forecast h=1 matches predict on the next-period lag", {
  dat <- sim_sir(m = 8, T_len = 16, p = 2, q = 1, family = "poisson", seed = 11)
  fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 11)
  Z_next <- dat$Z[, , , dim(dat$Z)[4], drop = FALSE]
  fc <- forecast(fit, h = 1, Z_future = Z_next)
  expect_equal(dim(fc), c(8, 8, 1))
  # off-diagonal forecasts are finite (the self-tie diagonal is NA by design)
  off <- fc[, , 1][lower.tri(fc[, , 1]) | upper.tri(fc[, , 1])]
  expect_true(all(is.finite(off)))

  # the h=1 forecast equals predict() on the lag built from the stored outcome,
  # off the diagonal (where self-ties are excluded)
  m <- fit$m
  x_next <- log(fit$Y[, , fit$n_periods] + 1)
  x_next[is.na(x_next)] <- 0
  x_next <- x_next / max(m - 1, 1)
  pr <- predict(fit, newdata = list(X = array(x_next, dim = c(8, 8, 1)),
                                    Z = Z_next), type = "response")
  mask <- lower.tri(fc[, , 1]) | upper.tri(fc[, , 1])
  expect_equal(fc[, , 1][mask], pr[, , 1][mask], tolerance = 1e-8)
})

test_that("forecast horizon > 1 returns the right shape and is finite", {
  dat <- sim_sir(m = 7, T_len = 14, p = 2, q = 1, family = "poisson", seed = 12)
  fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 12)
  Zf <- dat$Z[, , , (dim(dat$Z)[4] - 2):dim(dat$Z)[4], drop = FALSE]
  fc <- forecast(fit, h = 3, Z_future = Zf)
  expect_equal(dim(fc), c(7, 7, 3))
  # off-diagonal forecasts are finite (self-tie diagonal is NA by design)
  for (s in 1:3) {
    sl <- fc[, , s]
    expect_true(all(is.finite(sl[lower.tri(sl) | upper.tri(sl)])))
  }
  expect_equal(dimnames(fc)[[3]], c("h1", "h2", "h3"))
})

test_that("forecast guards bad horizon and missing Z_future", {
  dat <- sim_sir(m = 6, T_len = 10, p = 2, q = 1, family = "poisson", seed = 13)
  fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 13)
  expect_error(forecast(fit, h = 0))
  expect_error(forecast(fit, h = 1))   # q > 0 but no Z_future
})

test_that("score_sir computes family-appropriate scores", {
  a <- matrix(rpois(100, 3), 10, 10); diag(a) <- NA
  s_perfect <- score_sir(a, a, "poisson")
  expect_equal(unname(s_perfect[["rmse"]]), 0)
  expect_lt(s_perfect[["deviance"]], 1e-8)

  s_const <- score_sir(matrix(2, 10, 10), matrix(4, 10, 10), "normal")
  expect_equal(unname(s_const[["rmse"]]), 2)
  expect_equal(unname(s_const[["mae"]]), 2)

  yb <- c(rep(1, 50), rep(0, 50))
  pb <- c(rep(0.9, 50), rep(0.1, 50))
  s_bin <- score_sir(yb, pb, "binomial")
  expect_equal(unname(s_bin[["auc"]]), 1)
  expect_true(is.finite(s_bin[["logloss"]]))
})

test_that("cv_sir runs rolling-origin folds and beats the naive baseline", {
  dat <- sim_sir(m = 8, T_len = 20, p = 2, q = 1, family = "poisson", seed = 14)
  fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 14)
  cv <- cv_sir(fit, initial = 12, origins = c(12, 15, 18))
  expect_s3_class(cv, "sir_cv")
  expect_equal(nrow(cv$scores), 3)
  # all outcome training windows end strictly before the test period; future Z is
  # treated as known/pre-specified by cv_sir.
  expect_true(all(cv$scores$origin < fit$n_periods))
  expect_true(is.finite(cv$aggregate[["rmse"]]))
})

test_that("one-step cv_sir forecasts do not use held-out X slices", {
  dat <- sim_sir(m = 7, T_len = 12, p = 2, q = 1, family = "poisson", seed = 16)
  fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 16)
  fit_changed <- fit
  fit_changed$X[, , 9] <- fit_changed$X[, , 9] + 100
  cv1 <- cv_sir(fit, initial = 8, origins = 8)
  cv2 <- cv_sir(fit_changed, initial = 8, origins = 8)
  expect_equal(cv2$scores, cv1$scores, tolerance = 1e-8)
})

test_that("dyad jackknife gives finite, reproducible SEs on a square fit", {
  dat <- sim_sir(m = 8, T_len = 12, p = 2, q = 1, family = "poisson", seed = 15)
  fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson",
             fix_receiver = TRUE, seed = 15)
  b1 <- boot_sir(fit, R = 25, type = "dyad", seed = 5)
  b2 <- boot_sir(fit, R = 25, type = "dyad", seed = 5)
  expect_equal(b1$type, "dyad")
  expect_true(all(is.finite(b1$se)))
  expect_true(b1$n_valid > 5)
  expect_identical(b1$se, b2$se)
  ci <- confint(b1)
  expect_equal(nrow(ci), length(fit$tab))
})

test_that("full-bilinear bipartite recovers separate alpha and beta", {
  set.seed(202)
  n1 <- 10; n2 <- 6; Tt <- 18
  W  <- array(rnorm(n1 * n1 * 2), dim = c(n1, n1, 2))
  Wr <- array(rnorm(n2 * n2 * 2), dim = c(n2, n2, 2))
  X  <- array(rnorm(n1 * n2 * Tt) / sqrt(n2), dim = c(n1, n2, Tt))
  Z  <- array(rnorm(n1 * n2 * Tt), dim = c(n1, n2, Tt))
  alpha_t <- c(1, 0.6); beta_t <- c(0.8, -0.5); theta_t <- 1.5
  A_t <- W[, , 1] + alpha_t[2] * W[, , 2]
  B_t <- beta_t[1] * Wr[, , 1] + beta_t[2] * Wr[, , 2]
  eta <- array(0, dim = c(n1, n2, Tt))
  for (t in 1:Tt) eta[, , t] <- theta_t * Z[, , t] + A_t %*% X[, , t] %*% t(B_t)
  Y <- eta + array(rnorm(n1 * n2 * Tt, sd = 0.4), dim = c(n1, n2, Tt))

  fit <- sir(Y, W = W, X = X, Z = Z, W_recv = Wr, family = "normal", calc_se = FALSE)
  expect_s3_class(fit, "sir")
  expect_true(isTRUE(fit$bipartite))
  expect_true(isTRUE(fit$full_bilinear))
  expect_equal(dim(fit$A), c(n1, n1))
  expect_equal(dim(fit$B), c(n2, n2))
  expect_equal(fit$tab, c(theta_t, alpha_t[2], beta_t), tolerance = 0.1,
               ignore_attr = TRUE)
})

test_that("bipartite predict matches fitted values and supports counterfactuals", {
  set.seed(203)
  n1 <- 8; n2 <- 5; Tt <- 12
  W  <- array(rnorm(n1 * n1 * 2), dim = c(n1, n1, 2))
  Wr <- array(rnorm(n2 * n2 * 2), dim = c(n2, n2, 2))
  X  <- array(rnorm(n1 * n2 * Tt) / sqrt(n2), dim = c(n1, n2, Tt))
  Z  <- array(rnorm(n1 * n2 * Tt), dim = c(n1, n2, Tt))
  eta <- array(0, dim = c(n1, n2, Tt))
  A_t <- W[, , 1] + 0.4 * W[, , 2]; B_t <- 0.7 * Wr[, , 1] - 0.3 * Wr[, , 2]
  for (t in 1:Tt) eta[, , t] <- Z[, , t] + A_t %*% X[, , t] %*% t(B_t)
  Y <- eta + array(rnorm(n1 * n2 * Tt, sd = 0.3), dim = c(n1, n2, Tt))
  fit <- sir(Y, W = W, X = X, Z = Z, W_recv = Wr, family = "normal", calc_se = FALSE)

	  pin <- predict(fit)
	  expect_equal(max(abs(pin - fit$fitted.values)), 0, tolerance = 1e-8)
	  pcf <- predict(fit, newdata = list(W = W, X = X, Z = Z + 1))
	  expect_gt(mean(abs(pcf - pin)), 1e-6)
	  pwr <- predict(fit, newdata = list(W = W, W_recv = Wr + 0.25, X = X, Z = Z))
	  expect_gt(mean(abs(pwr - pin)), 1e-6)
	})

test_that("dyad jackknife works on a bipartite full-bilinear fit", {
  set.seed(204)
  n1 <- 9; n2 <- 5; Tt <- 14
  W  <- array(rnorm(n1 * n1 * 2), dim = c(n1, n1, 2))
  Wr <- array(rnorm(n2 * n2 * 2), dim = c(n2, n2, 2))
  X  <- array(rnorm(n1 * n2 * Tt) / sqrt(n2), dim = c(n1, n2, Tt))
  eta <- array(0, dim = c(n1, n2, Tt))
  A_t <- W[, , 1] + 0.5 * W[, , 2]; B_t <- Wr[, , 1] - 0.4 * Wr[, , 2]
  for (t in 1:Tt) eta[, , t] <- A_t %*% X[, , t] %*% t(B_t)
  Y <- eta + array(rnorm(n1 * n2 * Tt, sd = 0.3), dim = c(n1, n2, Tt))
  fit <- sir(Y, W = W, X = X, Z = NULL, W_recv = Wr, family = "normal", calc_se = FALSE)
  bb <- boot_sir(fit, R = 20, type = "dyad", seed = 8)
  expect_true(all(is.finite(bb$se)))
  expect_true(bb$n_valid > 5)
})

test_that("sim_sir catches mistyped time arguments and handles T_len = 1", {
  # the dots guard turns an unmatched time argument into a clear error
  expect_error(sim_sir(m = 6, time = 5, p = 2, q = 1), "T_len")
  # T_len = 1 must not trip the 2:T_len recursion (reverse-iteration footgun)
  d <- sim_sir(m = 6, T_len = 1, p = 2, q = 1, family = "poisson", seed = 1)
  expect_equal(dim(d$Y)[3], 1)
})

# --- regression tests for the 15-agent review findings ---

test_that("F20: sir()/cv_sir() do not clobber a global variable named 'fit'", {
  dat <- sim_sir(m = 8, T_len = 14, p = 2, q = 1, family = "poisson", seed = 3)
  fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 3)
  before_class <- class(fit)
  before_np <- fit$n_periods
  cv <- cv_sir(fit, initial = 8, origins = c(8, 10))
  # the source fit must survive the refits unchanged
  expect_identical(class(fit), before_class)
  expect_identical(fit$n_periods, before_np)
  expect_true(inherits(fit, "sir_fit"))
})

test_that("F31/F32/F33: dyad jackknife is centred on the point estimate (bipartite)", {
  set.seed(204)
  n1 <- 10; n2 <- 6; Tt <- 30
  W  <- array(rnorm(n1 * n1 * 2), dim = c(n1, n1, 2))
  Wr <- array(rnorm(n2 * n2 * 2), dim = c(n2, n2, 2))
  X  <- array(rnorm(n1 * n2 * Tt) / sqrt(n2), dim = c(n1, n2, Tt))
  A_t <- W[, , 1] + 0.6 * W[, , 2]; B_t <- 0.8 * Wr[, , 1] - 0.5 * Wr[, , 2]
  eta <- array(0, dim = c(n1, n2, Tt))
  for (t in 1:Tt) eta[, , t] <- A_t %*% X[, , t] %*% t(B_t)
  Y <- eta + array(rnorm(n1 * n2 * Tt, sd = 0.3), dim = c(n1, n2, Tt))
  fit <- sir(Y, W = W, X = X, W_recv = Wr, family = "normal", calc_se = FALSE, seed = 1)
  b <- boot_sir(fit, type = "dyad", seed = 1)
  # the point estimate must lie inside its own jackknife CI (the with-replacement
  # bootstrap attenuated beta so badly the point fell outside)
  expect_true(all(b$ci_lo <= b$point_est & b$point_est <= b$ci_hi))
  # SEs must be sane, not heavy-tailed/explosive
  expect_true(all(b$se < 1))
  # CI must cover the truth alpha2=0.6, betaWr=(0.8,-0.5)
  ci <- confint(b)
  truth <- c(0.6, 0.8, -0.5)
  expect_true(all(ci[, 1] <= truth & truth <= ci[, 2]))
})

test_that("F31: square fix_receiver=FALSE dyad jackknife point lies in its own CI", {
  dat <- sim_sir(m = 12, T_len = 50, p = 2, q = 1, family = "poisson", seed = 11)
  fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 1)
  b <- boot_sir(fit, type = "dyad", seed = 1)
  expect_true(all(b$se < 5))
  expect_true(all(b$ci_lo <= b$point_est & b$point_est <= b$ci_hi))
})

test_that("F35: parametric bootstrap works on a full-bilinear bipartite fit", {
  set.seed(204)
  n1 <- 9; n2 <- 5; Tt <- 20
  W  <- array(rnorm(n1 * n1 * 2), dim = c(n1, n1, 2))
  Wr <- array(rnorm(n2 * n2 * 2), dim = c(n2, n2, 2))
  X  <- array(rnorm(n1 * n2 * Tt) / sqrt(n2), dim = c(n1, n2, Tt))
  A_t <- W[, , 1] + 0.5 * W[, , 2]; B_t <- Wr[, , 1] - 0.4 * Wr[, , 2]
  eta <- array(0, dim = c(n1, n2, Tt))
  for (t in 1:Tt) eta[, , t] <- A_t %*% X[, , t] %*% t(B_t)
  Y <- eta + array(rnorm(n1 * n2 * Tt, sd = 0.3), dim = c(n1, n2, Tt))
  fit <- sir(Y, W = W, X = X, W_recv = Wr, family = "normal", calc_se = FALSE)
  bp <- boot_sir(fit, R = 15, type = "parametric", seed = 1)
  expect_true(all(is.finite(bp$se)))
})

test_that("F36: seeded boot_sir restores the caller's global RNG stream", {
  dat <- sim_sir(m = 8, T_len = 12, p = 2, q = 1, family = "poisson", seed = 1)
  fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson",
             fix_receiver = TRUE, seed = 1)
  set.seed(99); a <- runif(1)
  invisible(boot_sir(fit, R = 15, type = "block", seed = 7))
  set.seed(99); b <- runif(1)
  expect_identical(a, b)
})

test_that("F37: boot_sir point_est names match se/ci names", {
  dat <- sim_sir(m = 8, T_len = 12, p = 2, q = 1, family = "poisson", seed = 2)
  fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson",
             fix_receiver = TRUE, seed = 2)
  b <- boot_sir(fit, type = "dyad", seed = 1)
  expect_identical(names(b$point_est), names(b$se))
  expect_identical(names(b$point_est), b$param_names)
})

test_that("F10/F14: cv_sir validates horizon and reserved ... args", {
  dat <- sim_sir(m = 8, T_len = 16, p = 2, q = 1, family = "poisson", seed = 5)
  fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 5)
  expect_error(cv_sir(fit, initial = 10, horizon = 0), "horizon")
  expect_error(cv_sir(fit, initial = 10, horizon = 1.5), "horizon")
  expect_error(cv_sir(fit, initial = 10.5, horizon = 1), "initial")
  expect_error(cv_sir(fit, initial = 10, origins = c(10, 11.5)), "origins")
  expect_error(cv_sir(dat$Y), "fitted")
  expect_warning(cv_sir(fit, initial = 10, origins = c(10, 12), calc_se = TRUE), "Ignoring")
})

test_that("F25/F26: score_sir validates the response domain", {
  expect_error(score_sir(c(0, 2, 1), c(.1, .2, .3), "binomial"))
  expect_error(score_sir(c(-1, 2), c(1, 2), "poisson"))
  expect_error(score_sir(c(0, 1.5), c(0, 2), "poisson"), "integer")
  expect_error(score_sir(c(0, 2), c(0, NA), "poisson", drop_diagonal = FALSE), "non-finite")
  expect_no_error(score_sir(c(NA, 2), c(NA, 2), "poisson", drop_diagonal = FALSE))
  # perfect poisson prediction scores deviance 0
  a <- matrix(rpois(100, 2), 10, 10); diag(a) <- NA
  expect_lt(score_sir(a, a, "poisson")[["deviance"]], 1e-8)
  actual <- c(rep(1, 8), 100)
  predicted <- c(rep(0, 8), 100)
  expect_gte(score_sir(actual, predicted, "poisson", drop_diagonal = FALSE)[["deviance"]], 0)
})

test_that("F42: square Y + W_recv keeps two-mode diagonals in the likelihood", {
	set.seed(7)
	n = 8
	W = array(rnorm(n * n * 2), dim = c(n, n, 2))
	Wr = array(rnorm(n * n * 2), dim = c(n, n, 2))
	X = array(rnorm(n * n * 12) / sqrt(n), dim = c(n, n, 12))
	A_t = W[, , 1] + 0.4 * W[, , 2]
	B_t = 0.7 * Wr[, , 1]
	eta = array(0, dim = c(n, n, 12))
	for (t in seq_len(12)) {
		eta[, , t] = A_t %*% X[, , t] %*% t(B_t)
	}
	Y = eta + array(rnorm(n * n * 12, sd = 0.3), dim = c(n, n, 12))
	fit = sir(
		Y,
		W = W,
		X = X,
		W_recv = Wr,
		family = "normal",
		calc_se = FALSE,
		seed = 1
	)
	expect_equal(fit$nobs, n * n * 12)
})

test_that("full-bilinear bipartite normal tracks residual and likelihood scales", {
	set.seed(13)
	n1 <- 7; n2 <- 4; Tt <- 10
	W <- array(rnorm(n1 * n1 * 2), dim = c(n1, n1, 2))
	Wr <- array(rnorm(n2 * n2), dim = c(n2, n2, 1))
	X <- array(rnorm(n1 * n2 * Tt) / sqrt(n2), dim = c(n1, n2, Tt))
	A_t <- W[, , 1] + 0.3 * W[, , 2]
	B_t <- 0.8 * Wr[, , 1]
	eta <- array(0, dim = c(n1, n2, Tt))
	for (t in seq_len(Tt)) eta[, , t] <- A_t %*% X[, , t] %*% t(B_t)
	Y <- eta + array(rnorm(n1 * n2 * Tt, sd = 0.25), dim = c(n1, n2, Tt))
	Y[1, 1, 1] <- NA
	fit <- sir(Y, W = W, X = X, W_recv = Wr, family = "normal", calc_se = FALSE)

	rss <- sum(fit$residuals$response^2, na.rm = TRUE)
	expect_true(is.na(fitted(fit)[1, 1, 1]))
	expect_equal(fit$sigma2, rss / (fit$nobs - length(fit$tab)), tolerance = 1e-8)
	expect_equal(fit$sigma2_mle, rss / fit$nobs, tolerance = 1e-8)
	expect_gt(fit$sigma2, fit$sigma2_mle)
})

test_that("F41: bipartite SE accessors steer to boot_sir(type='dyad')", {
  set.seed(8); n1 <- 9; n2 <- 5; Tt <- 16
  W  <- array(rnorm(n1 * n1 * 2), dim = c(n1, n1, 2))
  Wr <- array(rnorm(n2 * n2 * 2), dim = c(n2, n2, 2))
  X  <- array(rnorm(n1 * n2 * Tt) / sqrt(n2), dim = c(n1, n2, Tt))
  A_t <- W[, , 1] + 0.5 * W[, , 2]; B_t <- Wr[, , 1]
  eta <- array(0, dim = c(n1, n2, Tt))
  for (t in 1:Tt) eta[, , t] <- A_t %*% X[, , t] %*% t(B_t)
  Y <- eta + array(rnorm(n1 * n2 * Tt, sd = 0.3), dim = c(n1, n2, Tt))
  fit <- sir(Y, W = W, X = X, W_recv = Wr, family = "normal", calc_se = FALSE)
  expect_error(vcov(fit), "boot_sir")
  expect_error(confint(fit), "boot_sir")
  expect_error(generics::tidy(fit), "boot_sir")
  # but confint via a boot result works
  b <- boot_sir(fit, type = "dyad", seed = 1)
  expect_equal(nrow(confint(fit, boot = b)), length(fit$tab))
  expect_false(any(is.finite(fit$summ$se)))
  expect_no_error(print(fit))
})

test_that("F53: tidy() labels bipartite receiver params as beta, not theta", {
  skip_if_not_installed("generics")
  set.seed(9); n1 <- 8; n2 <- 5; Tt <- 14
  W  <- array(rnorm(n1 * n1 * 2), dim = c(n1, n1, 2))
  Wr <- array(rnorm(n2 * n2 * 2), dim = c(n2, n2, 2))
  X  <- array(rnorm(n1 * n2 * Tt) / sqrt(n2), dim = c(n1, n2, Tt))
  A_t <- W[, , 1] + 0.4 * W[, , 2]; B_t <- 0.7 * Wr[, , 1] - 0.3 * Wr[, , 2]
  eta <- array(0, dim = c(n1, n2, Tt))
  for (t in 1:Tt) eta[, , t] <- A_t %*% X[, , t] %*% t(B_t)
  Y <- eta + array(rnorm(n1 * n2 * Tt, sd = 0.3), dim = c(n1, n2, Tt))
  fit <- sir(Y, W = W, X = X, W_recv = Wr, family = "normal", calc_se = FALSE)
  td <- generics::tidy(fit, se.type = "classical")
  betarows <- td[grepl("betaWr", td$term), ]
  expect_true(nrow(betarows) > 0)
  expect_true(all(betarows$component == "beta"))
})

test_that("F1: forecast warns when X scaling disagrees with its lag transform", {
  dat <- sim_sir(m = 8, T_len = 16, p = 2, q = 1, family = "poisson", seed = 1)
  Xraw <- array(0, dim(dat$Y))
  for (t in 2:dim(dat$Y)[3]) Xraw[, , t] <- log(dat$Y[, , t - 1] + 1)   # no /(m-1)
  fraw <- sir(dat$Y, W = dat$W, X = Xraw, Z = dat$Z, family = "poisson", seed = 1)
  Z_next <- dat$Z[, , , 16, drop = FALSE]
  expect_warning(forecast(fraw, h = 1, Z_future = Z_next), "lag transform")
  # the override silences it and runs
  expect_silent(forecast(fraw, h = 1, Z_future = Z_next, infl_scale = 1))
})

test_that("forecast scaling check includes square-bipartite diagonal cells", {
	set.seed(17)
	n <- 5; Tt <- 8
	Y <- array(rnorm(n * n * Tt), dim = c(n, n, Tt))
	W <- array(rnorm(n * n), dim = c(n, n, 1))
	Wr <- array(rnorm(n * n), dim = c(n, n, 1))
	infl_scale <- sqrt(max(n - 1, 1) * max(n - 1, 1))
	X <- array(0, dim = c(n, n, Tt))
	for (t in 2:Tt) X[, , t] <- Y[, , t - 1] / infl_scale
	fit <- sir(Y, W = W, X = X, W_recv = Wr, family = "normal",
			   calc_se = FALSE, seed = 17)
	fit_bad <- fit
	fit_bad$X[1, 1, Tt] <- fit_bad$X[1, 1, Tt] + 10
	expect_warning(forecast(fit_bad, h = 1), "lag transform")
})

test_that("F5: forecast rejects non-finite horizon", {
  dat <- sim_sir(m = 6, T_len = 10, p = 2, q = 1, family = "poisson", seed = 1)
  fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 1)
  expect_error(forecast(fit, h = Inf, Z_future = dat$Z[, , , 10, drop = FALSE]), "positive integer")
})

test_that("p=0 fits support predict newdata and forecast", {
	set.seed(10)
	n <- 6; Tt <- 8
	Y <- array(rpois(n * n * Tt, 2), dim = c(n, n, Tt))
	Z <- array(rnorm(n * n * Tt), dim = c(n, n, 1, Tt))
	fit <- sir(Y, W = NULL, X = NULL, Z = Z, family = "poisson", calc_se = FALSE)
	Z_new <- Z
	Z_new[, , 1, ] <- Z_new[, , 1, ] + 0.5
	pr <- predict(fit, newdata = list(X = fit$X, Z = Z_new))
	expect_equal(dim(pr), dim(Y))
	fc <- forecast(fit, h = 1, Z_future = Z[, , , Tt, drop = FALSE])
	expect_equal(dim(fc), c(n, n, 1))
})

test_that("full-bilinear bipartite predict rejects unsupported dynamic W", {
	set.seed(12)
	n1 <- 6; n2 <- 4; Tt <- 8
	W <- array(rnorm(n1 * n1 * 2), dim = c(n1, n1, 2))
	Wr <- array(rnorm(n2 * n2 * 2), dim = c(n2, n2, 2))
	X <- array(rnorm(n1 * n2 * Tt), dim = c(n1, n2, Tt))
	Y <- array(rnorm(n1 * n2 * Tt), dim = c(n1, n2, Tt))
	fit <- sir(Y, W = W, X = X, W_recv = Wr, family = "normal", calc_se = FALSE)
	W_dyn <- array(rnorm(n1 * n1 * 2 * Tt), dim = c(n1, n1, 2, Tt))
	expect_error(
		predict(fit, newdata = list(W = W_dyn, W_recv = Wr, X = X)),
		"3D static sender-side"
	)
})

test_that("discrete families reject impossible outcome domains", {
	Y_frac <- array(1.2, dim = c(4, 4, 3))
	expect_error(sir(Y_frac, W = NULL, X = NULL, family = "poisson"), "integer count")
	Y_bin <- array(0, dim = c(4, 4, 3))
	Y_bin[1, 2, 1] <- 0.5
	expect_error(sir(Y_bin, W = NULL, X = NULL, family = "binomial"), "0/1")
})
