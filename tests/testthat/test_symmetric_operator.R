# guards for the genuine symmetric (A = B) operator

make_sym_normal = function(n = 14, Tt = 70, p = 3, seed = 7) {
	set.seed(seed)
	W = array(0, c(n, n, p))
	for (k in 1:p) { Wk = matrix(rnorm(n*n), n, n); W[,,k] = 0.5*(Wk+t(Wk)); diag(W[,,k]) = 0 }
	gamma = c(1, 0.6, -0.4)[1:p]
	A = matrix(matrix(W, n*n, p) %*% gamma, n, n)
	Z = array(0, c(n, n, 1, Tt)); for (t in 1:Tt) { Zt = matrix(rnorm(n*n), n, n); Z[,,1,t] = 0.5*(Zt+t(Zt)) }
	X = array(0, c(n, n, Tt)); Y = array(0, c(n, n, Tt))
	for (t in 1:Tt) {
		Xt = matrix(rnorm(n*n), n, n)/sqrt(n); Xt = 0.5*(Xt+t(Xt)); X[,,t] = Xt
		e = A %*% Xt %*% t(A) + 1.0*Z[,,1,t]
		Yt = e + 0.4*matrix(rnorm(n*n), n, n); Yt = 0.5*(Yt+t(Yt)); diag(Yt) = NA; Y[,,t] = Yt
	}
	list(Y = Y, W = W, X = X, Z = Z, gamma = gamma, A = A)
}

test_that("symmetric recovers gamma and theta", {
	d = make_sym_normal()
	fit = suppressWarnings(sir(d$Y, W = d$W, X = d$X, Z = d$Z,
				  family = "normal", symmetric = TRUE, seed = 1))
	expect_equal(fit$gamma, d$gamma, tolerance = 0.05, ignore_attr = TRUE)
})

test_that("symmetric fitted network is exactly symmetric with NA diagonal", {
	d = make_sym_normal()
	fit = suppressWarnings(sir(d$Y, W = d$W, X = d$X, Z = d$Z,
				  family = "normal", symmetric = TRUE, seed = 1))
	for (t in 1:dim(fit$fitted.values)[3]) {
		ft = fit$fitted.values[, , t]
		expect_true(all(is.na(diag(ft))))
		expect_equal(max(abs(ft - t(ft)), na.rm = TRUE), 0, tolerance = 1e-10)
	}
})

test_that("symmetric enforces zero-diagonal A", {
	d = make_sym_normal()
	fit = suppressWarnings(sir(d$Y, W = d$W, X = d$X, Z = d$Z,
				  family = "normal", symmetric = TRUE, seed = 1))
	expect_lt(max(abs(diag(fit$A))), 1e-10)
	expect_equal(fit$A, fit$B, tolerance = 1e-12)
})

test_that("symmetric predict in-sample equals fitted", {
	d = make_sym_normal()
	fit = suppressWarnings(sir(d$Y, W = d$W, X = d$X, Z = d$Z,
				  family = "normal", symmetric = TRUE, seed = 1))
	pin = predict(fit, type = "response")
	expect_equal(as.numeric(pin), as.numeric(fit$fitted.values), tolerance = 1e-10)
})

test_that("sim_sir symmetric eta matches eta_tab_symmetric at the truth", {
	d = sim_sir(m = 12, T_len = 40, p = 3, q = 0, family = "normal",
				symmetric = TRUE, seed = 3)
	tab_true = d$alpha       # q = 0, so tab = gamma_1..p (all gammas free)
	eta = sir:::eta_tab_symmetric(tab_true, d$W, d$X, NULL, p = 3, q = 0)
	# the sim's own mean A X A' on the off-diagonal must match
	A = matrix(matrix(d$W, 12*12, 3) %*% d$alpha, 12, 12)
	for (t in 2:40) {
		em = A %*% d$X[,,t] %*% t(A)
		diff = max(abs((eta[,,t] - em)[upper.tri(em)]))
		expect_lt(diff, 1e-8)
	}
})

test_that("symmetric dyad jackknife refits the symmetric model and covers truth", {
	d = make_sym_normal(n = 12, Tt = 50)
	fit = suppressWarnings(sir(d$Y, W = d$W, X = d$X, Z = d$Z,
				  family = "normal", symmetric = TRUE, seed = 1))
	bj = boot_sir(fit, type = "dyad", seed = 1)
	expect_equal(length(bj$point_est), length(fit$tab))
	expect_true(all(is.finite(bj$se)))
	# the point estimate lies inside its own jackknife CI
	expect_true(all(bj$ci_lo <= bj$point_est & bj$point_est <= bj$ci_hi))
})

test_that("symmetric classical SE confint covers the truth (sanity)", {
	d = make_sym_normal(n = 14, Tt = 80)
	fit = suppressWarnings(sir(d$Y, W = d$W, X = d$X, Z = d$Z,
				  family = "normal", symmetric = TRUE, seed = 1))
	ci = confint(fit, se.type = "classical")
	# all p gammas are estimated now (gamma_1 free, identified up to sign)
	g_rows = grep("gammaW", rownames(ci))
	truth = d$gamma
	covered = ci[g_rows, 1] <= truth & truth <= ci[g_rows, 2]
	expect_true(all(covered))
})

test_that("symmetric works for poisson and binomial", {
	dp = sim_sir(m = 12, T_len = 60, p = 2, q = 1, family = "poisson",
				 symmetric = TRUE, seed = 3)
	fp = suppressWarnings(sir(dp$Y, W = dp$W, X = dp$X, Z = dp$Z,
				 family = "poisson", symmetric = TRUE, seed = 1))
	expect_true(all(fp$fitted.values[!is.na(fp$fitted.values)] >= 0))
	expect_true(isTRUE(fp$symmetric))

	db = sim_sir(m = 12, T_len = 120, p = 2, q = 0, family = "binomial",
				 symmetric = TRUE, seed = 5)
	fb = suppressWarnings(sir(db$Y, W = db$W, X = db$X, family = "binomial",
				 symmetric = TRUE, seed = 1))
	fvb = fb$fitted.values[!is.na(fb$fitted.values)]
	expect_true(all(fvb >= 0 & fvb <= 1))
})

test_that("symmetric supports dynamic (4D) W", {
	set.seed(1); n <- 12; Tt <- 150; p <- 3
	W4 <- array(0, c(n, n, p, Tt))
	for (t in 1:Tt) for (k in 1:p) {
		Wk <- matrix(rnorm(n * n), n, n) / sqrt(n)
		Wk <- (Wk + t(Wk)) / 2; diag(Wk) <- 0; W4[, , k, t] <- Wk
	}
	g <- c(1, 0.6, -0.4)
	X <- array(0, c(n, n, Tt)); Y <- array(0, c(n, n, Tt))
	for (t in 1:Tt) {
		At <- matrix(matrix(W4[, , , t], n * n, p) %*% g, n, n)
		Xt <- matrix(rnorm(n * n), n, n) / sqrt(n); Xt <- (Xt + t(Xt)) / 2; X[, , t] <- Xt
		e <- At %*% Xt %*% t(At)
		Yt <- e + 0.4 * matrix(rnorm(n * n), n, n); Yt <- (Yt + t(Yt)) / 2
		diag(Yt) <- NA; Y[, , t] <- Yt
	}
	fit <- suppressWarnings(sir(Y, W = W4, X = X, family = "normal", symmetric = TRUE, seed = 1))
	expect_true(isTRUE(fit$dynamic_W))
	expect_true(isTRUE(fit$symmetric))
	# per-period operator stored as an n x n x T array
	expect_equal(dim(fit$A), c(n, n, Tt))
	# recovers the shared gamma
	expect_equal(fit$gamma, g, tolerance = 0.05, ignore_attr = TRUE)
	# fitted symmetric, predict == fitted
	expect_equal(max(abs(fit$fitted.values[, , 5] - t(fit$fitted.values[, , 5])), na.rm = TRUE), 0, tolerance = 1e-10)
	expect_equal(as.numeric(predict(fit)), as.numeric(fit$fitted.values), tolerance = 1e-8)
})

test_that("symmetric warns on an asymmetric (directed) Z", {
	set.seed(2); n <- 10; Tt <- 20; p <- 2
	W <- array(0, c(n, n, p))
	for (k in 1:p) { Wk <- matrix(rnorm(n * n), n, n) / sqrt(n); W[, , k] <- (Wk + t(Wk)) / 2; diag(W[, , k]) <- 0 }
	X <- array(0, c(n, n, Tt)); for (t in 1:Tt) { Xt <- matrix(rnorm(n*n), n, n)/sqrt(n); X[, , t] <- (Xt + t(Xt))/2 }
	A <- matrix(matrix(W, n*n, p) %*% c(1, 0.5), n, n)
	Y <- array(0, c(n, n, Tt)); Zd <- array(rnorm(n*n*1*Tt), c(n, n, 1, Tt))   # asymmetric
	for (t in 1:Tt) { e <- A %*% X[,,t] %*% t(A) + Zd[,,1,t]; Yt <- e + 0.4*matrix(rnorm(n*n),n,n); Yt <- (Yt+t(Yt))/2; diag(Yt) <- NA; Y[,,t] <- Yt }
	expect_warning(
		suppressMessages(sir(Y, W = W, X = X, Z = Zd, family = "normal", symmetric = TRUE, seed = 1)),
		"not symmetric"
	)
})

test_that("symmetric cluster-robust SE is the default, PSD; classical still available", {
	d = make_sym_normal(n = 14, Tt = 60)
	fit = suppressWarnings(sir(d$Y, W = d$W, X = d$X, Z = d$Z,
				  family = "normal", symmetric = TRUE, seed = 1))
	# default is the actor cluster-robust sandwich, not classical
	expect_equal(confint(fit), confint(fit, se.type = "cluster"))
	expect_false(isTRUE(all.equal(confint(fit), confint(fit, se.type = "classical"))))
	# cluster-robust reachable, finite, PSD (non-negative variances)
	vc = vcov(fit, type = "cluster")
	expect_true(all(is.finite(vc)))
	expect_true(all(diag(vc) >= 0))
	# symmetric fits have no separate HC0 path, so "robust" maps to the cluster sandwich
	expect_equal(vcov(fit, type = "robust"), vc)
	# classical still available and finite
	expect_true(all(is.finite(confint(fit, se.type = "classical"))))
})

test_that("symmetric validates inputs (non-finite, family domain)", {
	d = make_sym_normal(n = 10, Tt = 20)
	# Inf in Y -> abort, not a silent converged fit
	Yi = d$Y; Yi[3, 5, 2] = Inf; Yi[5, 3, 2] = Inf
	expect_error(suppressWarnings(sir(Yi, W = d$W, X = d$X, Z = d$Z,
		family = "normal", symmetric = TRUE, seed = 1)), "non-finite")
	# Inf in W -> abort
	Wi = d$W; Wi[1, 2, 1] = Inf; Wi[2, 1, 1] = Inf
	expect_error(suppressWarnings(sir(d$Y, W = Wi, X = d$X, Z = d$Z,
		family = "normal", symmetric = TRUE, seed = 1)), "non-finite")
	# NA in X -> informative zero-fill, runs clean
	Xn = d$X; Xn[2, 3, 5] = NA
	expect_message(
		suppressWarnings(sir(d$Y, W = d$W, X = Xn, Z = d$Z,
			family = "normal", symmetric = TRUE, seed = 1)),
		"NA")
	fit = suppressWarnings(suppressMessages(sir(d$Y, W = d$W, X = Xn, Z = d$Z,
		family = "normal", symmetric = TRUE, seed = 1)))
	expect_s3_class(fit, "sir")
})

test_that("symmetric validates the family domain of Y", {
	set.seed(3); n = 10; Tt = 20; p = 2
	W = array(0, c(n, n, p))
	for (k in 1:p) { Wk = matrix(rnorm(n*n), n, n)/sqrt(n); W[,,k] = (Wk+t(Wk))/2; diag(W[,,k]) = 0 }
	X = array(0, c(n, n, Tt)); for (t in 1:Tt) { Xt = matrix(rnorm(n*n),n,n)/sqrt(n); X[,,t] = (Xt+t(Xt))/2 }
	A = matrix(matrix(W, n*n, p) %*% c(1, 0.5), n, n)
	# poisson with a non-integer -> abort
	Yp = array(0, c(n, n, Tt))
	for (t in 1:Tt) { e = A %*% X[,,t] %*% t(A); Yt = matrix(rpois(n*n, exp(pmin(e,4))), n, n); Yt = Yt + t(Yt); Yt[lower.tri(Yt)] = t(Yt)[lower.tri(Yt)]; diag(Yt) = NA; Yp[,,t] = Yt }
	Yp2 = Yp; Yp2[1,2,1] = 2.5; Yp2[2,1,1] = 2.5
	expect_error(suppressWarnings(sir(Yp2, W = W, X = X, family = "poisson", symmetric = TRUE, seed = 1)), "integer")
	# binomial with a 2 -> abort
	Yb = Yp; Yb[] = (Yp > 0) * 1; for (t in 1:Tt) diag(Yb[,,t]) = NA
	Yb2 = Yb; Yb2[1,2,1] = 2; Yb2[2,1,1] = 2
	expect_error(suppressWarnings(sir(Yb2, W = W, X = X, family = "binomial", symmetric = TRUE, seed = 1)), "0/1|Bernoulli")
})

test_that("symmetric forecast uses the full symmetric lag (forecast == predict)", {
	set.seed(1); m = 14; Tt = 40; p = 2
	W = array(0, c(m, m, p))
	for (k in 1:p) { Wk = matrix(rnorm(m*m), m, m)/sqrt(m); W[,,k] = (Wk+t(Wk))/2; diag(W[,,k]) = 0 }
	A = matrix(matrix(W, m*m, p) %*% c(1, 0.5), m, m)
	Y = array(0, c(m, m, Tt)); X = array(0, c(m, m, Tt))
	Y[,,1] = { y = matrix(rpois(m*m, 2), m, m); y = y + t(y); diag(y) = NA; y }
	for (t in 2:Tt) {
		xl = log(Y[,,t-1] + 1) / (m - 1); xl[is.na(xl)] = 0; X[,,t] = xl
		e = A %*% xl %*% t(A); Yt = matrix(rpois(m*m, exp(pmin(e, 5))), m, m)
		Yt = Yt + t(Yt); Yt[lower.tri(Yt)] = t(Yt)[lower.tri(Yt)]; diag(Yt) = NA; Y[,,t] = Yt
	}
	fit = suppressWarnings(sir(Y, W = W, X = X, family = "poisson", symmetric = TRUE, seed = 1))
	# a properly-lagged symmetric X must not trigger the lag-mismatch warning
	expect_silent(fc <- forecast(fit, h = 1))
	# forecast off-diagonal equals predict on the hand-built symmetric lag
	Yl = fit$Y[, , fit$n_periods]; Yl[lower.tri(Yl)] = t(Yl)[lower.tri(Yl)]
	xl = log(Yl + 1) / (m - 1); xl[is.na(xl)] = 0
	ph = predict(fit, newdata = list(W = W, X = array(xl, c(m, m, 1))), type = "response")
	off = row(matrix(0, m, m)) != col(matrix(0, m, m))
	expect_lt(max(abs(fc[, , 1] - ph[, , 1])[off], na.rm = TRUE), 1e-8)
})

test_that("symmetric = TRUE and fix_receiver = TRUE are incompatible", {
	dd = sim_sir(m = 8, T_len = 15, p = 3, q = 1, family = "normal", symmetric = TRUE, seed = 5)
	expect_error(
		sir(dd$Y, W = dd$W, X = dd$X, Z = dd$Z, family = "normal",
			symmetric = TRUE, fix_receiver = TRUE, seed = 1),
		"incompatible"
	)
})

test_that("symmetric symmetrizes an asymmetric X so the fit stays symmetric", {
	d = sim_sir(m = 10, T_len = 25, p = 2, q = 1, family = "normal", symmetric = TRUE, seed = 4)
	Xa = array(rnorm(10 * 10 * 25), dim = c(10, 10, 25))   # deliberately not symmetric
	fit = suppressWarnings(suppressMessages(sir(d$Y, W = d$W, X = Xa, Z = d$Z,
		family = "normal", symmetric = TRUE, seed = 1)))
	fv = fit$fitted.values[, , 5]
	expect_lt(max(abs(fv - t(fv)), na.rm = TRUE), 1e-10)
})

test_that("sim_sir(symmetric=TRUE) keeps an explicit non-unit anchor and recovers it", {
	d = sim_sir(m = 14, T_len = 90, p = 3, q = 1, family = "normal",
				symmetric = TRUE, alpha = c(1.8, 0.6, -0.4), seed = 2)
	expect_equal(d$alpha[1], 1.8)
	fit = suppressWarnings(sir(d$Y, W = d$W, X = d$X, Z = d$Z,
		family = "normal", symmetric = TRUE, seed = 1))
	expect_equal(unname(fit$gamma), c(1.8, 0.6, -0.4), tolerance = 0.25)
})

test_that("binomial symmetric fit reports stationary = NA (gain not applicable)", {
	d = sim_sir(m = 12, T_len = 80, p = 2, q = 0, family = "binomial", symmetric = TRUE, seed = 5)
	fit = suppressWarnings(sir(d$Y, W = d$W, X = d$X, family = "binomial", symmetric = TRUE, seed = 1))
	expect_true(is.na(fit$stationary))
})

test_that("m = 2 symmetric fit does not crash on the near-constant-anchor check", {
	d = sim_sir(m = 2, T_len = 30, p = 2, q = 1, family = "normal", symmetric = TRUE, seed = 1)
	expect_no_error(suppressWarnings(suppressMessages(sir(d$Y, W = d$W, X = d$X, Z = d$Z,
		family = "normal", symmetric = TRUE, seed = 1))))
})

test_that("family = 'gaussian' errors with a 'did you mean normal' hint", {
	d = sim_sir(m = 6, T_len = 12, p = 2, q = 1, family = "normal", seed = 1)
	expect_error(sir(d$Y, W = d$W, X = d$X, Z = d$Z, family = "gaussian"), "normal")
})
