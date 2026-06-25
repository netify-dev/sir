# standard-error paths, bootstrap, and symmetric-operator edge cases

make_sym = function(m = 8, T_len = 22, fam = "normal", seed = 5) {
	dat = sim_sir(m = m, T_len = T_len, p = 3, q = 1, family = fam,
				  symmetric = TRUE, seed = seed)
	suppressWarnings(sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z,
						 family = fam, symmetric = TRUE, seed = 1))
}

test_that("boot_sir block and parametric work for symmetric fits", {
	fs = make_sym()
	bb = suppressWarnings(boot_sir(fs, R = 8, type = "block", seed = 1, trace = FALSE))
	bp = suppressWarnings(boot_sir(fs, R = 8, type = "parametric", seed = 1, trace = FALSE))
	expect_gte(bb$n_valid, 2)
	expect_gte(bp$n_valid, 2)
	expect_true(all(is.finite(confint(fs, boot = bb))))
	expect_true(all(is.finite(confint(fs, boot = bp))))
})

test_that("parametric bootstrap still works for directed fits (regression)", {
	dat = sim_sir(m = 8, T_len = 22, p = 3, q = 1, family = "normal", seed = 5)
	fd = suppressWarnings(sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z,
							  family = "normal", seed = 1))
	bp = suppressWarnings(boot_sir(fd, R = 8, type = "parametric", seed = 1, trace = FALSE))
	expect_gte(bp$n_valid, 2)
})

test_that("dynamic (4D) W: default SE accessors use the cluster-robust sandwich", {
	m = 8; T_len = 15; p = 3
	Wt = array(0, dim = c(m, m, p, T_len))
	set.seed(11)
	for (t in seq_len(T_len)) for (k in seq_len(p)) {
		Wk = matrix(rnorm(m * m), m, m) / sqrt(m)
		Wk = (Wk + t(Wk)) / 2; diag(Wk) = 0
		Wt[, , k, t] = Wk
	}
	dd = sim_sir(m = m, T_len = T_len, p = 3, q = 1, family = "normal",
				 symmetric = TRUE, seed = 9)
	fdyn = suppressWarnings(sir(dd$Y, W = Wt, X = dd$X, Z = dd$Z,
							    family = "normal", symmetric = TRUE, seed = 1))
	# analytic SEs are available, so the jackknife fallback did not fire
	expect_null(fdyn$se_source)
	# default no-arg accessors all return finite SEs from the cluster sandwich
	suppressMessages({
		expect_true(all(is.finite(vcov(fdyn))))
		expect_true(all(is.finite(confint(fdyn))))
		expect_true(all(is.finite(tidy(fdyn, conf.int = TRUE)$conf.low)))
	})
	# the default IS the cluster-robust sandwich (not classical) for dynamic W
	suppressMessages(v_def <- vcov(fdyn))
	expect_equal(v_def, vcov(fdyn, type = "cluster"))
	expect_equal(attr(vcov(fdyn, type = "cluster"), "cluster_df"), m - 1L)
	# an explicit cluster request now succeeds instead of aborting
	suppressMessages(Vc <- vcov(fdyn, type = "cluster"))
	expect_true(all(is.finite(Vc)))
	# classical is still reachable and differs from the cluster sandwich
	se_cl = sqrt(diag(vcov(fdyn, type = "classical")))
	se_cr = sqrt(diag(Vc))
	expect_false(isTRUE(all.equal(se_cl, se_cr)))
})

test_that("tidy() cluster p-values use the same t(G-1) reference as confint()", {
	fs = make_sym()
	td = tidy(fs)
	G = attr(vcov(fs), "cluster_df") + 1L
	se = sqrt(diag(vcov(fs)))
	tstat = coef(fs) / se
	p_t = 2 * pt(abs(tstat), df = G - 1, lower.tail = FALSE)
	expect_equal(unname(td$p.value), unname(p_t), tolerance = 1e-10)
	# tidy conf bounds match confint() exactly (same SE + same reference)
	ci = confint(fs)
	expect_equal(unname(cbind(td2 <- tidy(fs, conf.int = TRUE))$conf.low), unname(ci[, 1]), tolerance = 1e-10)
	expect_equal(unname(td2$conf.high), unname(ci[, 2]), tolerance = 1e-10)
})

test_that("directed tidy() cluster p-values also use t(G-1)", {
	dat = sim_sir(m = 8, T_len = 22, p = 3, q = 1, family = "normal", seed = 5)
	fd = suppressWarnings(sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z,
							  family = "normal", seed = 1))
	td = tidy(fd)
	G = attr(vcov(fd), "cluster_df") + 1L
	se = sqrt(diag(vcov(fd)))
	p_t = 2 * pt(abs(coef(fd) / se), df = G - 1, lower.tail = FALSE)
	expect_equal(unname(td$p.value), unname(p_t), tolerance = 1e-10)
})

test_that("symmetric recovers relative gammas when the anchor loading != 1", {
	m = 14; T = 110
	W = array(0, dim = c(m, m, 3)); set.seed(99)
	for (k in 1:3) { Wk = matrix(rnorm(m * m), m, m) / sqrt(m); W[, , k] = (Wk + t(Wk)) / 2; diag(W[, , k]) = 0 }
	Z = array(0, dim = c(m, m, 1, T)); for (t in 1:T) { Zt = matrix(rnorm(m * m), m, m); Z[, , 1, t] = (Zt + t(Zt)) / 2 }
	gt = c(2, 1.2, -0.8)                      # anchor loading 2, NOT 1
	A = matrix(matrix(W, m * m, 3) %*% gt, m, m)
	X = array(0, dim = c(m, m, T)); Y = array(0, dim = c(m, m, T)); set.seed(7)
	for (t in 1:T) {
		if (t == 1) Xt = matrix(0, m, m) else { Xt = Y[, , t - 1]; Xt[is.na(Xt)] = 0; Xt = (Xt + t(Xt)) / 2; Xt = Xt / (m - 1) }
		X[, , t] = Xt; eta = A %*% Xt %*% t(A) + 1.2 * Z[, , 1, t]
		Yt = eta + matrix(rnorm(m * m, sd = 0.4), m, m); Yt = (Yt + t(Yt)) / 2; diag(Yt) = NA; Y[, , t] = Yt
	}
	fit = suppressWarnings(sir(Y, W = W, X = X, Z = Z, family = "normal", symmetric = TRUE, seed = 1))
	# all p gammas estimated, recovering the true (non-unit) anchor and loadings
	expect_equal(length(fit$gamma), 3L)
	expect_equal(unname(fit$gamma), gt, tolerance = 0.15)
	# sign convention: gamma_1 >= 0
	expect_gte(fit$gamma[1], 0)
})

test_that("bipartite cluster SE counts senders and receivers as distinct actors", {
	set.seed(1); n1 = 10; n2 = 7; Tb = 12
	Yb = array(rpois(n1 * n2 * Tb, 2), dim = c(n1, n2, Tb))
	Wb = array(rnorm(n1 * n1 * 3), dim = c(n1, n1, 3)); for (k in 1:3) diag(Wb[, , k]) = 0
	Xb = array(0, dim = c(n1, n2, Tb)); for (t in 2:Tb) Xb[, , t] = log(Yb[, , t - 1] + 1)
	fb = suppressWarnings(sir(Yb, W = Wb, X = Xb, family = "poisson",
							  fix_receiver = TRUE, bipartite = TRUE, seed = 1))
	expect_equal(attr(vcov(fb), "cluster_df"), n1 + n2 - 1L)
	# one-mode square clusters node i once (sender == receiver)
	dd = sim_sir(m = 12, T_len = 25, p = 3, q = 1, family = "normal", seed = 5)
	fd = suppressWarnings(sir(dd$Y, W = dd$W, X = dd$X, Z = dd$Z, family = "normal", seed = 1))
	expect_equal(attr(vcov(fd), "cluster_df"), 12L - 1L)
})

test_that("predict() errors on bare W/X/Z instead of silently ignoring them", {
	dd = sim_sir(m = 10, T_len = 20, p = 3, q = 1, family = "normal", symmetric = TRUE, seed = 5)
	fit = suppressWarnings(sir(dd$Y, W = dd$W, X = dd$X, Z = dd$Z, family = "normal", symmetric = TRUE, seed = 1))
	expect_error(predict(fit, W = dd$W), "cannot be passed directly")
	# the correct newdata path still works
	expect_silent(p1 <- predict(fit, newdata = list(W = dd$W, X = dd$X, Z = dd$Z)))
})

test_that("symmetric direct effects carry the (Z) tag and gamma rows cover all p", {
	dd = sim_sir(m = 10, T_len = 20, p = 3, q = 1, family = "normal", symmetric = TRUE, seed = 5)
	fit = suppressWarnings(sir(dd$Y, W = dd$W, X = dd$X, Z = dd$Z, family = "normal", symmetric = TRUE, seed = 1))
	rn = rownames(fit$summ)
	expect_true(any(grepl("^\\(Z\\)", rn)))
	expect_equal(sum(grepl("\\(gammaW\\)", rn)), 3L)
	expect_true(isTRUE(fit$stationary))      # a well-behaved fit is stationary
})

test_that("sir auto-falls-back to jackknife SEs when analytic SEs are unavailable", {
	# full-bilinear bipartite has no closed-form covariance; sir() should attach a
	# delete-one-actor jackknife covariance automatically so the accessors work.
	set.seed(204)
	n1 = 10; n2 = 6; Tn = 40
	Wf = array(rnorm(n1 * n1 * 2), c(n1, n1, 2))
	Wr = array(rnorm(n2 * n2 * 2), c(n2, n2, 2))
	Xf = array(rnorm(n1 * n2 * Tn) / sqrt(n2), c(n1, n2, Tn))
	At = Wf[, , 1] + 0.6 * Wf[, , 2]
	Bt = 0.8 * Wr[, , 1] - 0.5 * Wr[, , 2]
	Yf = array(0, c(n1, n2, Tn))
	for (t in 1:Tn) Yf[, , t] = At %*% Xf[, , t] %*% t(Bt) +
		matrix(rnorm(n1 * n2, sd = 0.3), n1, n2)

	fit = suppressMessages(sir(Yf, W = Wf, X = Xf, W_recv = Wr,
							   family = "normal", seed = 1))

	expect_identical(fit$se_source, "jackknife")
	# vcov() / confint() now return values instead of aborting
	V = vcov(fit)
	expect_true(is.matrix(V) && all(is.finite(diag(V))))
	ci = confint(fit)
	expect_true(is.matrix(ci) && nrow(ci) == length(coef(fit)))
	# the jackknife intervals cover the generating weights
	truth = c(0.6, 0.8, -0.5)
	expect_true(all(ci[, 1] <= truth & truth <= ci[, 2]))

	# summary() carries se_source so it prints the jackknife (not classical) footnote
	expect_identical(summary(fit)$se_source, "jackknife")

	# the attached covariance matches boot_sir(type="dyad")'s own (summed-margin)
	# covariance, so vcov(fit) agrees with the explicit resampling path
	bs = suppressWarnings(boot_sir(fit, type = "dyad"))
	expect_equal(unname(diag(V)), unname(diag(bs$cov)), tolerance = 1e-2)

	# tidy() returns finite SEs on the fallback fit for EVERY se.type -- the
	# robust path used to return all-NA std.error and a false "refit with
	# calc_se = TRUE" message
	for (st in c("cluster", "classical", "robust")) {
		td = suppressMessages(tidy(fit, se.type = st))
		expect_true(all(is.finite(td$std.error)),
					info = paste("tidy se.type =", st))
	}

	# a clean directed fit keeps its analytic SEs (no fallback)
	dd = sim_sir(m = 14, T_len = 70, p = 2, q = 1, family = "normal", seed = 7)
	f2 = suppressMessages(sir(dd$Y, W = dd$W, X = dd$X, Z = dd$Z,
							  family = "normal", seed = 1))
	expect_null(f2$se_source)
	expect_true(all(is.finite(sqrt(diag(vcov(f2))))))
})
