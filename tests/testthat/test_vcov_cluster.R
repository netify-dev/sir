# check cluster-robust variance calculations

test_that("per-observation scores reconstruct the kernel score outer-product", {
	set.seed(1)
	dat = sim_sir(m = 10, T_len = 20, p = 2, q = 2, family = "poisson", seed = 1)
	fit = sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z,
			  family = "poisson", calc_se = TRUE, seed = 1)

	sc = sir:::.sir_obs_scores(fit)
	Z_list = sir:::prepare_Z_list(dat$Z)
	gH = sir:::cpp_mll_gH(fit$tab, dat$Y, dat$W, dat$X, Z_list, "poisson")

	# crossprod(scores) == kernel shess (the BHHH meat); colSums(scores) == -grad.
	expect_equal(dim(crossprod(sc$scores)), dim(gH$shess))
	expect_lt(max(abs(crossprod(sc$scores) - gH$shess)) / max(abs(gH$shess)), 1e-8)
	expect_lt(max(abs(colSums(sc$scores) + as.numeric(gH$grad))), 1e-5)
})

test_that("cluster-robust vcov is a valid PSD covariance, exposed via vcov/confint", {
	set.seed(2)
	dat = sim_sir(m = 10, T_len = 20, p = 2, q = 2, family = "poisson", seed = 2)
	fit = sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z,
			  family = "poisson", calc_se = TRUE, seed = 2)

	V = vcov(fit, type = "cluster")
	expect_equal(dim(V), dim(fit$vcov))
	expect_true(all(is.finite(V)))
	expect_gt(min(eigen(V, symmetric = TRUE, only.values = TRUE)$values), -1e-8)
	expect_true(isSymmetric(unname(V), tol = 1e-8))

	ci = confint(fit, se.type = "cluster")
	expect_equal(nrow(ci), length(coef(fit)))
	expect_true(all(is.finite(ci)))

	# dyad clustering also produces a valid covariance
	Vd = vcov(fit, type = "dyad")
	expect_true(all(is.finite(Vd)))
})

test_that("cluster and twoway aliases agree", {
	set.seed(3)
	dat = sim_sir(m = 9, T_len = 18, p = 2, q = 2, family = "poisson", seed = 3)
	fit = sir(
		dat$Y,
		W = dat$W,
		X = dat$X,
		Z = dat$Z,
		family = "poisson",
		calc_se = TRUE,
		seed = 3
	)

	expect_equal(vcov(fit, type = "cluster"), vcov(fit, type = "twoway"))
	expect_true(all(is.finite(vcov(fit, type = "dyad"))))
})
