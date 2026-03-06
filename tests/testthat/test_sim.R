test_that("sim_sir generates valid Poisson network data", {
	dat = sim_sir(m = 8, T_len = 5, p = 2, q = 1,
				 family = "poisson", seed = 42)

	expect_equal(dim(dat$Y), c(8, 8, 5))
	expect_equal(dim(dat$W), c(8, 8, 2))
	expect_equal(dim(dat$X), c(8, 8, 5))
	expect_equal(dim(dat$Z), c(8, 8, 1, 5))
	expect_equal(length(dat$alpha), 2)
	expect_equal(dat$alpha[1], 1)
	expect_equal(length(dat$beta), 2)
	expect_equal(length(dat$theta), 1)
	expect_true(all(dat$Y >= 0))
	# diagonal should be 0
	for (t in 1:5) expect_true(all(diag(dat$Y[,,t]) == 0))
})

test_that("sim_sir generates valid binomial data", {
	dat = sim_sir(m = 6, T_len = 3, p = 2, q = 0,
				 family = "binomial", seed = 42)

	expect_null(dat$Z)
	expect_true(all(dat$Y %in% c(0, 1)))
})

test_that("sim_sir generates valid normal data", {
	dat = sim_sir(m = 6, T_len = 3, p = 2, q = 2,
				 family = "normal", sigma = 2, seed = 42)

	expect_equal(dim(dat$Z), c(6, 6, 2, 3))
	# normal data can be negative
	expect_true(any(dat$Y < 0))
})

test_that("sim_sir respects user-supplied parameters", {
	dat = sim_sir(m = 5, T_len = 3, p = 2, q = 1,
				 family = "poisson",
				 alpha = c(1, 0.5), beta = c(0.3, -0.2),
				 theta = c(0.1), seed = 42)

	expect_equal(dat$alpha, c(1, 0.5))
	expect_equal(dat$beta, c(0.3, -0.2))
	expect_equal(dat$theta, 0.1)
})

test_that("sim_sir output can be fit by sir()", {
	dat = sim_sir(m = 8, T_len = 10, p = 2, q = 1,
				 family = "poisson", seed = 42)

	fit = sir(dat$Y, dat$W, dat$X, dat$Z,
			 family = "poisson", fix_receiver = TRUE,
			 calc_se = FALSE, max_iter = 5)

	expect_s3_class(fit, "sir")
	expect_equal(length(coef(fit)), 1 + 2)  # q=1 theta + p=2 alpha
})
