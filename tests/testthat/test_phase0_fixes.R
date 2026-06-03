# tests for seed handling, convergence reporting, and robust standard errors

test_that("seed argument produces reproducible results", {
	dat = sim_sir(m = 10, T_len = 12, p = 2, q = 1, family = "poisson", seed = 123)
	fit1 = sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 7)
	fit2 = sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 7)
	expect_equal(fit1$tab, fit2$tab)
	expect_equal(fit1$summ$se, fit2$summ$se)
})

test_that("fix_receiver reports real convergence", {
	dat = sim_sir(m = 10, T_len = 12, p = 2, q = 1, family = "poisson", seed = 42)
	fit = sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", fix_receiver = TRUE, seed = 1)
	expect_true(is.logical(fit$convergence))
	expect_true(fit$convergence)
})

test_that("optim method reports convergence status", {
	dat = sim_sir(m = 10, T_len = 12, p = 2, q = 1, family = "poisson", seed = 11)
	fit = sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", method = "optim", seed = 1)
	expect_true(is.logical(fit$convergence))
})

test_that("normal logLik includes sigma^2 in df", {
	dat = sim_sir(m = 10, T_len = 12, p = 2, q = 1, family = "normal", seed = 5)
	fit = sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "normal", seed = 1)
	ll = logLik(fit)
	expect_true(is.finite(as.numeric(ll)))
	expect_true(attr(ll, "df") >= length(fit$tab) + 1)
})

test_that("fix_receiver robust SE is a non-degenerate sandwich", {
	dat = sim_sir(m = 12, T_len = 15, p = 2, q = 1, family = "poisson", seed = 1)
	fit_fr = sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", fix_receiver = TRUE, calc_se = TRUE, seed = 1)
	expect_true(all(is.finite(fit_fr$summ$rse)))
	expect_true(all(fit_fr$summ$rse > 0))
	expect_false(all(fit_fr$summ$rse == fit_fr$summ$se))
})

test_that("classical and robust SEs are both available and differ", {
	dat = sim_sir(m = 12, T_len = 15, p = 2, q = 1, family = "poisson", seed = 2)
	fit = sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", calc_se = TRUE, seed = 2)
	expect_true(all(is.finite(fit$summ$se)))
	expect_true(all(is.finite(fit$summ$rse)))
	expect_true(all(fit$summ$se > 0))
	expect_true(all(fit$summ$rse > 0))
	expect_false(isTRUE(all.equal(fit$summ$se, fit$summ$rse)))
})

test_that("seed restoration leaves global RNG untouched", {
	set.seed(999)
	before = .Random.seed
	dat = sim_sir(m = 8, T_len = 10, p = 2, q = 1, family = "poisson", seed = 5)
	fit = sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 5)
	after = .Random.seed
	expect_equal(before, after)
})

test_that("glm_converged helper handles speedglm and glm objects", {
	skip_if_not_installed("speedglm")
	dat = sim_sir(m = 10, T_len = 12, p = 2, q = 1, family = "poisson", seed = 3)
	fit = sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", fix_receiver = TRUE, seed = 1)
	expect_true(fit$convergence)
})
