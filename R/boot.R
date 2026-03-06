
#' Bootstrap Inference for SIR Model Parameters
#'
#' Computes bootstrap standard errors and confidence intervals for SIR model
#' parameters. This is the recommended approach for inference when the Hessian
#' is singular or ill-conditioned, which is common in models with bilinear
#' influence terms (i.e., when \code{fix_receiver = FALSE}).
#'
#' Two bootstrap strategies are available:
#' \describe{
#'   \item{block}{Resamples time periods with replacement. Preserves the
#'     within-period dependence structure. Best when T is moderately large
#'     (T >= 10).}
#'   \item{parametric}{Simulates new outcome arrays from the fitted model
#'     using the estimated parameters and the specified family distribution.
#'     Better when T is small but the model is well-specified.}
#' }
#'
#' Each replicate refits the full SIR model. Replicates that fail to converge
#' are dropped and reported. Standard errors are column standard deviations
#' of the successful replicates. Confidence intervals use the percentile
#' method.
#'
#' @param sir_fit A fitted \code{sir} object from \code{\link{sir}}.
#' @param R Integer. Number of bootstrap replicates. Default is 200. Increase
#'   to 500-1000 for publication-quality intervals.
#' @param type Character. Bootstrap type: \code{"block"} (default) resamples
#'   time periods with replacement; \code{"parametric"} simulates new outcomes
#'   from the fitted model.
#' @param seed Optional integer for reproducibility. Sets the random seed
#'   before resampling.
#' @param trace Logical. If TRUE, prints progress every 10 replicates.
#'
#' @return An object of class \code{"boot_sir"} with components:
#' \describe{
#'   \item{coefs}{R x n_params matrix of bootstrap coefficient estimates.
#'     Rows for failed replicates contain NA.}
#'   \item{se}{Named numeric vector of bootstrap standard errors (one per
#'     parameter).}
#'   \item{ci_lo}{Lower 2.5\% percentile bounds.}
#'   \item{ci_hi}{Upper 97.5\% percentile bounds.}
#'   \item{point_est}{Point estimates from the original fit.}
#'   \item{param_names}{Character vector of parameter names.}
#'   \item{n_valid}{Number of successful bootstrap replicates.}
#'   \item{n_total}{Total number of replicates attempted.}
#'   \item{type}{The bootstrap type used.}
#'   \item{family}{The distribution family.}
#' }
#'
#' @examples
#' \dontrun{
#' model <- sir(Y, W, X, family = "poisson")
#'
#' # block bootstrap with 200 replicates
#' boot_result <- boot_sir(model, R = 200, seed = 42)
#' print(boot_result)
#'
#' # use bootstrap CIs with confint
#' confint(model, boot = boot_result)
#'
#' # parametric bootstrap
#' boot_par <- boot_sir(model, R = 200, type = "parametric")
#' }
#'
#' @seealso \code{\link{confint.sir}} to use bootstrap intervals,
#'   \code{\link{confint.boot_sir}} for direct interval extraction.
#' @export
boot_sir <- function(sir_fit, R = 200, type = c("block", "parametric"),
					 seed = NULL, trace = FALSE) {
	type <- match.arg(type)
	if (!is.null(seed)) set.seed(seed)

	Y <- sir_fit$Y
	W <- sir_fit$W
	X <- sir_fit$X
	Z <- sir_fit$Z
	family <- sir_fit$family
	T_len <- dim(Y)[3]
	n_params <- length(sir_fit$tab)
	p <- if (is.null(W)) 0L else dim(W)[3]
	q <- if (is.null(Z)) 0L else dim(Z)[3]

	boot_coefs <- matrix(NA, R, n_params)
	colnames(boot_coefs) <- rownames(sir_fit$summ)

	for (b in 1:R) {
		if (trace && b %% 10 == 0) cli::cli_inform("Bootstrap {.val {b}}/{.val {R}}")

		if (type == "block") {
			# resample time periods with replacement
			t_idx <- sample(1:T_len, T_len, replace = TRUE)
			Y_b <- Y[,, t_idx, drop = FALSE]
			X_b <- X[,, t_idx, drop = FALSE]
			Z_b <- if (!is.null(Z) && length(dim(Z)) == 4) {
				Z[,,, t_idx, drop = FALSE]
			} else {
				Z
			}
		} else {
			# parametric: simulate from fitted model
			eta <- eta_tab(sir_fit$tab, W, X, Z)
			if (family == "normal") {
				sigma <- sqrt(sir_fit$sigma2)
				Y_b <- array(rnorm(length(eta), mean = eta, sd = sigma),
							 dim = dim(Y))
			} else if (family == "poisson") {
				lambda <- exp(eta)
				lambda[lambda > 1e6] <- 1e6
				Y_b <- array(rpois(length(eta), lambda = lambda),
							 dim = dim(Y))
			} else {
				prob <- 1 / (1 + exp(-eta))
				Y_b <- array(rbinom(length(eta), 1, prob),
							 dim = dim(Y))
			}
			# preserve diagonal NA structure
			for (tt in 1:T_len) diag(Y_b[,,tt]) <- NA
			X_b <- X
			Z_b <- Z
		}

		tryCatch({
			fit_b <- sir(Y_b, W, X_b, Z_b, family = family,
						 method = "ALS", calc_se = FALSE,
						 fix_receiver = isTRUE(sir_fit$fix_receiver),
						 kron_mode = isTRUE(sir_fit$kron_mode),
						 max_iter = 100, tol = 1e-6)
			tab_b <- fit_b$tab

			# sign alignment for bilinear models: bootstrap replicates
			# can converge to the reflected solution (-alpha, -beta).
			# flip influence params if they point away from the original.
			if (!isTRUE(sir_fit$fix_receiver) && p > 1 && length(tab_b) > q) {
				infl_idx <- (q + 1):length(tab_b)
				if (sum(sir_fit$tab[infl_idx] * tab_b[infl_idx]) < 0) {
					tab_b[infl_idx] <- -tab_b[infl_idx]
				}
			}

			boot_coefs[b, ] <- tab_b
		}, error = function(e) {
			# leave as NA for failed replicates
		})
	}

	# compute statistics from successful replicates
	valid <- apply(boot_coefs, 1, function(r) !any(is.na(r)))
	n_valid <- sum(valid)

	if (n_valid < 10) {
		cli::cli_warn("Only {.val {n_valid}} valid bootstrap replicates (of {.val {R}}). Results unreliable.")
	}

	boot_se <- apply(boot_coefs[valid, , drop = FALSE], 2, sd)
	boot_ci <- apply(boot_coefs[valid, , drop = FALSE], 2,
					 quantile, probs = c(0.025, 0.975))

	result <- list(
		coefs = boot_coefs,
		se = boot_se,
		ci_lo = boot_ci[1, ],
		ci_hi = boot_ci[2, ],
		point_est = sir_fit$tab,
		param_names = rownames(sir_fit$summ),
		n_valid = n_valid,
		n_total = R,
		type = type,
		family = sir_fit$family
	)
	class(result) <- "boot_sir"
	result
}

#' Print Bootstrap SIR Results
#'
#' Displays a table of point estimates, bootstrap standard errors, and 95\%
#' percentile confidence intervals.
#'
#' @param x A \code{boot_sir} object from \code{\link{boot_sir}}.
#' @param digits Number of significant digits. Default uses
#'   \code{getOption("digits") - 3}.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns the \code{boot_sir} object.
#' @export
print.boot_sir <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
	cli::cli_text("\n")
	cli::cli_text("{.strong Bootstrap SIR Results}")
	cli::cli_text("Type: {.field {x$type}} | Replicates: {.val {x$n_valid}}/{.val {x$n_total}} valid")

	# build results table
	tab <- cbind(
		Estimate = x$point_est,
		`Boot SE` = x$se,
		`2.5 %` = x$ci_lo,
		`97.5 %` = x$ci_hi
	)
	rownames(tab) <- x$param_names

	cli::cli_text("")
	print(round(tab, digits))
	cli::cli_text("")
	invisible(x)
}

#' Summary of Bootstrap SIR Results
#'
#' Displays detailed bootstrap results including the coefficient table,
#' significance indicators (whether the 95\% CI excludes zero), and
#' bootstrap distribution summaries (mean, sd, median).
#'
#' @param object A \code{boot_sir} object from \code{\link{boot_sir}}.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns the \code{boot_sir} object.
#' @export
summary.boot_sir <- function(object, ...) {
	cli::cli_h1("Bootstrap SIR Results")
	cli::cli_text("Type: {.field {object$type}} | Family: {.field {object$family}}")
	cli::cli_text("Replicates: {.val {object$n_valid}} valid of {.val {object$n_total}} total ({.val {round(100 * object$n_valid / object$n_total, 1)}}%)")

	cli::cli_rule()
	cli::cli_h2("Parameter Estimates")

	tab <- cbind(
		Estimate = object$point_est,
		`Boot SE` = object$se,
		`2.5 %` = object$ci_lo,
		`97.5 %` = object$ci_hi
	)
	rownames(tab) <- object$param_names

	# flag parameters where CI includes zero
	covers_zero <- object$ci_lo <= 0 & object$ci_hi >= 0
	sig <- ifelse(covers_zero, " ", "*")
	tab_print <- cbind(tab, ` ` = sig)

	print(noquote(format(round(tab, 4), width = 10)))
	cli::cli_text("{.emph * = 95% CI excludes zero}")

	# bootstrap distribution summary
	valid <- apply(object$coefs, 1, function(r) !any(is.na(r)))
	boot_valid <- object$coefs[valid, , drop = FALSE]

	cli::cli_rule()
	cli::cli_h2("Bootstrap Distribution")
	dist_tab <- data.frame(
		mean = colMeans(boot_valid),
		sd = apply(boot_valid, 2, sd),
		median = apply(boot_valid, 2, median),
		row.names = object$param_names
	)
	print(round(dist_tab, 4))

	invisible(object)
}

#' Confidence Intervals from Bootstrap SIR Results
#'
#' Computes percentile confidence intervals at the specified level from the
#' bootstrap coefficient distribution.
#'
#' @param object A \code{boot_sir} object from \code{\link{boot_sir}}.
#' @param parm Character vector of parameter names or integer indices. If
#'   NULL (default), returns intervals for all parameters.
#' @param level Confidence level between 0 and 1. Default is 0.95.
#' @param ... Additional arguments (unused).
#' @return A matrix with one row per parameter and columns for the lower
#'   and upper bounds, labeled by percentage.
#' @export
confint.boot_sir <- function(object, parm = NULL, level = 0.95, ...) {
	a <- (1 - level) / 2
	pct <- paste0(format(100 * c(a, 1 - a), trim = TRUE, digits = 3), " %")

	valid <- apply(object$coefs, 1, function(r) !any(is.na(r)))
	boot_valid <- object$coefs[valid, , drop = FALSE]
	ci <- t(apply(boot_valid, 2, quantile, probs = c(a, 1 - a)))
	rownames(ci) <- object$param_names
	colnames(ci) <- pct

	if (!is.null(parm)) {
		ci <- ci[parm, , drop = FALSE]
	}
	ci
}
