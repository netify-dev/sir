
# broom tidiers: tidy / glance / augment should produce well-formed frames that
# slot into modelsummary / gtsummary pipelines.

make_fit = function(seed = 1, fix_receiver = FALSE) {
	dat = sim_sir(m = 9, T_len = 8, p = 2, q = 2, family = "poisson", seed = seed)
	sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson",
		fix_receiver = fix_receiver, calc_se = TRUE, seed = seed)
}

test_that("tidy.sir returns the broom-standard columns and one row per parameter", {
	fit = make_fit()
	td = tidy(fit)
	expect_s3_class(td, "data.frame")
	expect_true(all(c("term", "component", "estimate", "std.error",
					  "statistic", "p.value") %in% names(td)))
	expect_equal(nrow(td), length(fit$tab))
	expect_true(all(is.finite(td$estimate)))
	# component tagging covers the three parameter blocks
	expect_true(all(td$component %in% c("theta", "alpha", "beta")))
	expect_true(all(c("theta", "alpha", "beta") %in% td$component))
})

test_that("tidy.sir conf.int adds interval columns and respects se.type", {
	fit = make_fit()
	td = tidy(fit, conf.int = TRUE)
	expect_true(all(c("conf.low", "conf.high") %in% names(td)))
	expect_true(all(td$conf.low <= td$conf.high))

	# robust SEs give a different std.error column than classical
	tdc = tidy(fit, se.type = "classical")
	tdr = tidy(fit, se.type = "robust")
	expect_false(isTRUE(all.equal(tdc$std.error, tdr$std.error)))

	td_cluster = tidy(fit, se.type = "cluster", conf.int = TRUE)
	expect_equal(nrow(td_cluster), length(fit$tab))
	expect_true(all(is.finite(td_cluster$std.error)))
	expect_true(all(td_cluster$conf.low <= td_cluster$conf.high))
})

test_that("glance.sir returns a one-row model summary", {
	fit = make_fit()
	gl = glance(fit)
	expect_s3_class(gl, "data.frame")
	expect_equal(nrow(gl), 1L)
	expect_true(all(c("nobs", "df", "logLik", "AIC", "BIC",
					  "family", "method", "converged") %in% names(gl)))
	expect_equal(gl$nobs, fit$nobs)
	expect_equal(gl$AIC, AIC(fit))
})

test_that("augment.sir returns long fitted/residual data with no likelihood-excluded cells", {
	fit = make_fit()
	au = augment(fit)
	expect_s3_class(au, "data.frame")
	expect_true(all(c("sender", "receiver", "time", ".observed",
					  ".fitted", ".resid") %in% names(au)))
	# diagonal/self-tie and masked cells must be dropped
	expect_true(all(!is.na(au$.observed)))
	# response residual identity holds
	expect_equal(au$.resid, au$.observed - au$.fitted, tolerance = 1e-8)
})

test_that("tidiers work for a fix_receiver model (no beta block)", {
	fit = make_fit(fix_receiver = TRUE)
	td = tidy(fit)
	expect_false("beta" %in% td$component)
	expect_true(all(c("theta", "alpha") %in% td$component))
	expect_equal(nrow(glance(fit)), 1L)
})
