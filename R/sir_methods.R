#' @importFrom cli cli_h1 cli_h2 cli_text cli_ul cli_alert_success cli_alert_warning cli_rule cli_abort cli_warn cli_inform
#' @importFrom stats printCoefmat deviance quantile pnorm symnum AIC BIC logLik fitted residuals coef confint nobs qnorm sd
#' @keywords internal
NULL

.sir_deparse_one <- function(x) paste(deparse(x), collapse = " ")

#' Extract Model Coefficients from a SIR Model
#'
#' Returns the estimated parameter values from a fitted SIR model. These
#' include exogenous covariate effects (theta), sender influence weights
#' (alpha), and receiver influence weights (beta).
#'
#' @param object A fitted \code{sir} object from \code{\link{sir}}.
#' @param ... Additional arguments (unused).
#' @return Named numeric vector of estimated coefficients. Names use the
#'   \code{(Z)} / \code{(alphaW)} / \code{(betaW)} tagging from
#'   \code{summary()}: \code{(Z)} = exogenous direct effects (theta);
#'   \code{(alphaW)} = sender-influence weights (alpha_2..alpha_p, since
#'   alpha_1 = 1 is fixed for identifiability); \code{(betaW)} = receiver-
#'   influence weights (beta).
#' @seealso \code{\link{confint.sir}} for confidence intervals,
#'   \code{\link{vcov.sir}} for the variance-covariance matrix.
#' @export
coef.sir <- function(object, ...) {
	out <- object$summ$coef
	# name the coefficients so coef() matches summary()/tidy() and supports
	# name-based subsetting
	if (!is.null(object$summ) && !is.null(rownames(object$summ))) {
	names(out) <- rownames(object$summ)
	}
	out
}

#' Extract Fitted Values from a SIR Model
#'
#' Returns the fitted values on the response scale: expected counts for Poisson,
#' probabilities for binomial, or conditional means for normal.
#'
#' @param object A fitted \code{sir} object from \code{\link{sir}}.
#' @param ... Additional arguments (unused).
#' @return An array with the same dimensions as \code{Y} (n1 x n2 x T)
#'   containing fitted values on the response scale.
#' @export
fitted.sir <- function(object, ...) {
	object$fitted.values
}

#' Extract Residuals from a SIR Model
#'
#' Returns residuals of the specified type. Response residuals are raw
#' (Y - fitted). Pearson residuals are standardized by the variance function.
#' Deviance residuals are signed square roots of the individual deviance
#' contributions and are most useful for diagnostic plots.
#'
#' @param object A fitted \code{sir} object from \code{\link{sir}}.
#' @param type Character string specifying residual type: \code{"deviance"}
#'   (default), \code{"pearson"}, or \code{"response"}.
#' @param ... Additional arguments (unused).
#' @return An array with the same dimensions as \code{Y} containing the
#'   requested residuals. Contains NA where \code{Y} is missing.
#' @export
residuals.sir <- function(object, type = c("deviance", "pearson", "response"), ...) {
	type <- match.arg(type)
	if (is.null(object$residuals)) return(NULL)
	object$residuals[[type]]
}

#' Extract Log-Likelihood from a SIR Model
#'
#' Returns the log-likelihood at convergence as a \code{logLik} object, with
#' attributes for degrees of freedom and number of observations. This allows
#' \code{AIC()} and \code{BIC()} to work directly on the result.
#'
#' @param object A fitted \code{sir} object from \code{\link{sir}}.
#' @param ... Additional arguments (unused).
#' @return A \code{logLik} object with attributes \code{df} (number of
#'   estimated parameters) and \code{nobs} (number of observations).
#' @export
logLik.sir <- function(object, ...) {
	val <- object$ll
	# the normal family estimates an extra scale parameter (sigma^2), which must
	# count toward df so AIC/BIC align with lm()/glm() conventions
	df <- length(object$summ$coef) + isTRUE(object$family == "normal")
	attr(val, "df") <- df
	attr(val, "nobs") <- object$nobs
	class(val) <- "logLik"
	val
}

#' Extract Number of Observations from a SIR Model
#'
#' Returns the number of non-missing dyad-time observations used in fitting.
#' For one-mode square networks, diagonal entries (self-loops) are excluded from
#' this count. Square bipartite fits keep diagonal cells because rows and columns
#' are distinct actor sets.
#'
#' @param object A fitted \code{sir} object from \code{\link{sir}}.
#' @param ... Additional arguments (unused).
#' @return Integer count of observations.
#' @export
nobs.sir <- function(object, ...) {
	object$nobs
}

#' Akaike Information Criterion for a SIR Model
#'
#' Computes AIC = -2 * log-likelihood + k * (number of parameters).
#' Use this to compare SIR models with different specifications (e.g.,
#' different numbers of influence covariates).
#'
#' @param object A fitted \code{sir} object from \code{\link{sir}}.
#' @param ... Additional arguments for comparison with other models.
#' @param k Numeric penalty per parameter (default 2 for standard AIC).
#' @return Numeric AIC value. Lower is better.
#' @export
AIC.sir <- function(object, ..., k = 2) {
	objects <- list(object, ...)
	if (length(objects) > 1L) {
		out <- data.frame(
			df = vapply(objects, function(x) attr(logLik(x), "df"), numeric(1)),
			AIC = vapply(objects, function(x) {
				ll <- logLik(x)
				-2 * as.numeric(ll) + k * attr(ll, "df")
			}, numeric(1))
		)
		names <- vapply(substitute(list(object, ...))[-1L], .sir_deparse_one, character(1))
		row.names(out) <- names
		return(out)
	}
	ll <- logLik(object)
	-2 * as.numeric(ll) + k * attr(ll, "df")
}

#' Bayesian Information Criterion for a SIR Model
#'
#' Computes BIC = -2 * log-likelihood + log(nobs) * (number of parameters).
#' BIC penalizes model complexity more heavily than AIC for large samples.
#'
#' @param object A fitted \code{sir} object from \code{\link{sir}}.
#' @param ... Additional arguments for comparison with other models.
#' @return Numeric BIC value. Lower is better.
#' @export
BIC.sir <- function(object, ...) {
	objects <- list(object, ...)
	if (length(objects) > 1L) {
		out <- data.frame(
			df = vapply(objects, function(x) attr(logLik(x), "df"), numeric(1)),
			BIC = vapply(objects, function(x) {
				ll <- logLik(x)
				nobs <- attr(ll, "nobs")
				-2 * as.numeric(ll) + log(nobs) * attr(ll, "df")
			}, numeric(1))
		)
		names <- vapply(substitute(list(object, ...))[-1L], .sir_deparse_one, character(1))
		row.names(out) <- names
		return(out)
	}
	ll <- logLik(object)
	nobs <- attr(ll, "nobs")
	-2 * as.numeric(ll) + log(nobs) * attr(ll, "df")
}

#' Summary of a Fitted SIR Model
#'
#' Produces a detailed summary of the fitted model including coefficient
#' estimates with standard errors and p-values, model fit statistics
#' (log-likelihood, AIC, BIC), convergence status, and summaries of the
#' estimated influence matrices A and B.
#'
#' @param object A fitted \code{sir} object from \code{\link{sir}}.
#' @param ... Additional arguments (unused).
#' @return An object of class \code{"summary.sir"} containing:
#'   \describe{
#'     \item{coefficients}{Data frame with columns \code{coef}, \code{se},
#'       \code{p.value}, and significance codes.}
#'     \item{loglik}{Log-likelihood at convergence.}
#'     \item{aic}{AIC value.}
#'     \item{bic}{BIC value.}
#'     \item{converged}{Logical convergence indicator.}
#'     \item{iterations}{Iteration count.}
#'     \item{A.summary}{List with mean, sd, and range of off-diagonal
#'       entries in the sender effects matrix.}
#'     \item{B.summary}{Same for the receiver effects matrix.}
#'   }
#' @seealso \code{\link{print.summary.sir}} for the printed output.
#' @export
summary.sir <- function(object, ...) {

	# create summary object
	ans <- list()
	ans$call <- object$call
	ans$family <- object$family
	ans$method <- object$method
	ans$symmetric <- isTRUE(object$symmetric)
	ans$fix_receiver <- isTRUE(object$fix_receiver)
	ans$se_reliable <- object$se_reliable
	ans$convergence <- isTRUE(object$convergence)

	# model dimensions
	ans$m <- if (!is.null(object$m)) object$m else nrow(object$A)
	ans$n1 <- if (!is.null(object$n1)) object$n1 else nrow(object$A)
	ans$n2 <- if (!is.null(object$n2)) object$n2 else ncol(object$B)
	ans$bipartite <- isTRUE(object$bipartite)
	ans$n_periods <- object$n_periods
	ans$p <- object$p
	ans$q <- object$q
	ans$nobs <- object$nobs

		# coefficients table with significance
		ans$coefficients <- object$summ
		if ("se" %in% colnames(ans$coefficients) &&
			isTRUE(object$convergence) &&
			!isFALSE(object$se_reliable) &&
			any(is.finite(ans$coefficients$se))) {
		z_scores <- ans$coefficients$coef / ans$coefficients$se
		ans$coefficients$p.value <- 2 * (1 - pnorm(abs(z_scores)))

	ans$coefficients$sig <- symnum(ans$coefficients$p.value,
									corr = FALSE, na = FALSE,
									cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
									symbols = c("***", "**", "*", ".", " "))
	}

	# flag full-bilinear bipartite so print steers inference to boot_sir(dyad)
	ans$full_bilinear <- isTRUE(object$full_bilinear)

	# model fit statistics
	ans$loglik <- object$ll
	ans$aic <- AIC(object)
	ans$bic <- BIC(object)
	ans$deviance <- if (!is.null(object$deviance)) object$deviance else -2 * object$ll
	ans$null.deviance <- object$null.deviance

	# residual deviance and (poisson/binomial) dispersion from observed vs fitted.
	# computed on the same cells the likelihood uses: drop NA and, for square
	# networks, the diagonal (self-loops are excluded from fitting). robust to
	# fix_receiver / symmetric / bipartite since it only reads Y and fitted.values.
	ans$resid.deviance <- NULL
	ans$dispersion <- NULL
	Y <- object$Y
	mu <- object$fitted.values
	if (!is.null(Y) && !is.null(mu) && all(dim(Y) == dim(mu))) {
	mask <- !is.na(Y) & !is.na(mu)
	# exclude the diagonal for square (directed/symmetric) networks
	if (!isTRUE(object$bipartite) && nrow(Y) == ncol(Y) && length(dim(Y)) == 3) {
	  diag_idx <- which(slice.index(Y, 1) == slice.index(Y, 2))
	  mask[diag_idx] <- FALSE
	}
	yv <- Y[mask]
	muv <- mu[mask]
	resid.df <- length(yv) - length(object$summ$coef)
	if (object$family == "poisson") {
	  dev_contrib <- 2 * (ifelse(yv > 0, yv * log(yv / muv), 0) - (yv - muv))
	  ans$resid.deviance <- sum(dev_contrib)
	  if (resid.df > 0) {
		ans$dispersion <- sum((yv - muv)^2 / muv) / resid.df
	  }
	} else if (object$family == "binomial") {
	  muc <- pmin(pmax(muv, 1e-10), 1 - 1e-10)
	  dev_contrib <- 2 * (ifelse(yv > 0, yv * log(yv / muc), 0) +
						  ifelse(yv < 1, (1 - yv) * log((1 - yv) / (1 - muc)), 0))
	  ans$resid.deviance <- sum(dev_contrib)
	  if (resid.df > 0) {
		ans$dispersion <- sum((yv - muc)^2 / (muc * (1 - muc))) / resid.df
	  }
	} else {
	  # normal: residual deviance is the residual sum of squares
	  ans$resid.deviance <- sum((yv - muv)^2)
	}
	}

	# convergence info
	ans$converged <- object$convergence
	ans$iterations <- object$iterations

	# influence matrices summary
	# handle dynamic W case where A/B are 3D arrays
	if (length(dim(object$A)) == 3) {
	# dynamic W: summarize across all time slices
	A_vals <- c()
	B_vals <- c()
	for (t in seq_len(dim(object$A)[3])) {
	  A_t <- object$A[,,t]
	  B_t <- object$B[,,t]
	  if (nrow(A_t) == ncol(A_t)) {
		A_vals <- c(A_vals, A_t[row(A_t) != col(A_t)])
	  } else {
		A_vals <- c(A_vals, c(A_t))
	  }
	  if (nrow(B_t) == ncol(B_t)) {
		B_vals <- c(B_vals, B_t[row(B_t) != col(B_t)])
	  } else {
		B_vals <- c(B_vals, c(B_t))
	  }
	}
	A_offdiag <- A_vals
	B_offdiag <- B_vals
	} else if (nrow(object$A) == ncol(object$A)) {
	A_offdiag <- object$A[row(object$A) != col(object$A)]
	B_offdiag <- if (nrow(object$B) == ncol(object$B)) {
	  object$B[row(object$B) != col(object$B)]
	} else {
	  c(object$B)
	}
	} else {
	A_offdiag <- c(object$A)
	B_offdiag <- c(object$B)
	}

	ans$A.summary <- list(
	mean = mean(A_offdiag),
	sd = sd(A_offdiag),
	range = range(A_offdiag)
	)

	ans$B.summary <- list(
	mean = mean(B_offdiag),
	sd = sd(B_offdiag),
	range = range(B_offdiag)
	)

	# sigma for normal family
	if (object$family == "normal" && !is.null(object$sigma2)) {
	ans$sigma <- sqrt(object$sigma2)
	}

	class(ans) <- "summary.sir"
	ans
}

#' Print a SIR Model Summary
#'
#' Displays the full model summary including network dimensions, family,
#' method, a coefficient table with significance stars (when standard errors
#' are available), model fit statistics, convergence status, and influence
#' matrix summaries.
#'
#' @param x A \code{summary.sir} object from \code{\link{summary.sir}}.
#' @param digits Number of significant digits to print. Default uses
#'   \code{getOption("digits") - 3}.
#' @param signif.stars Logical, whether to show significance stars beside
#'   p-values. Default uses \code{getOption("show.signif.stars")}.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns the summary object.
#' @export
print.summary.sir <- function(x, digits = max(3L, getOption("digits") - 3L),
							   signif.stars = getOption("show.signif.stars"), ...) {

	cli::cli_h1("Social Influence Regression Model")

	# network and model info
	if (isTRUE(x$bipartite)) {
	net_type <- "bipartite"
	cli::cli_text("{.strong Network:} {.val {x$n1}} senders x {.val {x$n2}} receivers, {.val {x$n_periods}} time periods ({net_type})")
	} else {
	net_type <- if (isTRUE(x$symmetric)) "symmetric (undirected)" else "directed"
	cli::cli_text("{.strong Network:} {.val {x$m}} nodes, {.val {x$n_periods}} time periods ({net_type})")
	}
	cli::cli_text("{.strong Family:} {.val {x$family}} | {.strong Method:} {.val {x$method}}")
	if (isTRUE(x$fix_receiver) && !isTRUE(x$symmetric)) {
	cli::cli_text("{.strong Receiver:} fixed (B = I)")
	}
	if (!is.null(x$nobs)) {
	cli::cli_text("{.strong Observations:} {.val {x$nobs}}")
	}

	cli::cli_rule()

	# coefficients
	cli::cli_h2("Coefficients")

	# treat an all-NA se column as "no SEs": full-bilinear bipartite fits carry an
	# se column of NAs, and steering users to vcov(type="robust") there dead-ends.
	have_se <- "se" %in% colnames(x$coefficients) &&
		"p.value" %in% colnames(x$coefficients) &&
		any(is.finite(x$coefficients$se))
	if (nrow(x$coefficients) > 0) {
	if (have_se) {
	  coef_mat <- as.matrix(x$coefficients[, c("coef", "se", "t_se", "p.value")])
	  colnames(coef_mat) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

	  if (signif.stars && "sig" %in% colnames(x$coefficients)) {
		printCoefmat(coef_mat, digits = digits, signif.stars = TRUE,
					 P.values = TRUE, has.Pvalue = TRUE)
	  } else {
		print(round(coef_mat, digits))
	  }
		  cli::cli_text("{.emph Std. errors above are classical (Hessian-based) and assume independent dyad-periods. For reported intervals, use {.code confint(fit)}; cluster-robust intervals are the default for supported static fits.}")
		} else {
	  # no SEs available
	  coef_mat <- as.matrix(x$coefficients[, "coef", drop = FALSE])
	  colnames(coef_mat) <- "Estimate"
	  print(round(coef_mat, digits))
		  if (isTRUE(x$full_bilinear)) {
			cli::cli_text("{.emph Analytic SEs are unavailable for full-bilinear bipartite fits. Use {.code boot_sir(fit, type = \"dyad\")} for inference.}")
		  } else if (isFALSE(x$se_reliable)) {
			cli::cli_text("{.emph Hessian-based SEs were marked unreliable for this fit, so p-values and stars are suppressed. Use {.code boot_sir()} or simplify/refit the model before reporting inference.}")
		  } else {
			cli::cli_text("{.emph (Standard errors not computed. Use calc_se = TRUE or boot_sir() for inference.)}")
		  }
	}
	} else {
	cli::cli_text("{.emph No coefficients estimated.}")
	}

	cli::cli_rule()

	# model fit
	cli::cli_h2("Model Fit")

	fit_stats <- c(
	paste0("Log-Likelihood: ", sprintf("%.2f", x$loglik)),
	paste0("AIC: ", sprintf("%.2f", x$aic)),
	paste0("BIC: ", sprintf("%.2f", x$bic))
	)

	if (!is.null(x$resid.deviance) && is.finite(x$resid.deviance)) {
	fit_stats <- c(fit_stats,
				   paste0("Residual deviance: ", sprintf("%.2f", x$resid.deviance)))
	}

	# dispersion > 1 flags overdispersion for poisson/binomial
	if (!is.null(x$dispersion) && is.finite(x$dispersion)) {
	fit_stats <- c(fit_stats,
				   paste0("Dispersion (Pearson chi-sq / resid df): ",
						  sprintf("%.3f", x$dispersion)))
	}

	if (!is.null(x$sigma)) {
	fit_stats <- c(fit_stats,
				   paste0("Residual std. error: ", sprintf("%.4f", x$sigma)))
	}

	cli::cli_ul(fit_stats)

	if (!is.null(x$dispersion) && is.finite(x$dispersion) && x$dispersion > 1.5) {
	cli::cli_alert_warning("Dispersion {sprintf('%.2f', x$dispersion)} > 1 suggests overdispersion relative to the {x$family} variance assumption.")
	}

	# convergence
	if (x$converged) {
	cli::cli_alert_success("Converged in {.val {x$iterations}} iterations")
	} else {
	cli::cli_alert_warning("Did not converge")
	}

	cli::cli_rule()

	# influence matrices summary
	if (x$p > 0) {
	if (isTRUE(x$symmetric)) {
	  cli::cli_h2("Influence Matrix")
	  cli::cli_text("{.strong A matrix (influence weights):}")
	} else {
	  cli::cli_h2("Influence Matrices")
	  cli::cli_text("{.strong A matrix (sender effects):}")
	}

	cli::cli_ul(c(
	  paste0("Mean: ", sprintf("%.4f", x$A.summary$mean)),
	  paste0("SD: ", sprintf("%.4f", x$A.summary$sd)),
	  paste0("Range: [", sprintf("%.4f", x$A.summary$range[1]), ", ",
			 sprintf("%.4f", x$A.summary$range[2]), "]")
	))

	if (!isTRUE(x$fix_receiver)) {
	  cli::cli_text("{.strong B matrix (receiver effects):}")
	  cli::cli_ul(c(
		paste0("Mean: ", sprintf("%.4f", x$B.summary$mean)),
		paste0("SD: ", sprintf("%.4f", x$B.summary$sd)),
		paste0("Range: [", sprintf("%.4f", x$B.summary$range[1]), ", ",
			   sprintf("%.4f", x$B.summary$range[2]), "]")
	  ))
	}
	}

	invisible(x)
}


#' Print a Fitted SIR Model
#'
#' Displays a compact overview of the fitted model: network dimensions,
#' family, method, convergence status, log-likelihood, AIC, and coefficient
#' estimates. Use \code{summary()} for a more detailed report with p-values
#' and influence matrix summaries.
#'
#' @param x A fitted \code{sir} object from \code{\link{sir}}.
#' @param digits Number of digits to print. Default uses
#'   \code{getOption("digits") - 3}.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns the sir object.
#' @export
print.sir <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
	cli::cli_text("\n")
	cli::cli_text("{.strong Social Influence Regression Model}")

	# network info
	if (isTRUE(x$bipartite)) {
	cli::cli_text("{.val {x$n1}} senders x {.val {x$n2}} receivers, {.val {x$n_periods}} time periods (bipartite)")
	} else {
	m <- if (!is.null(x$m)) x$m else nrow(x$A)
	net_type <- if (isTRUE(x$symmetric)) "symmetric (undirected)" else "directed"
	cli::cli_text("{.val {m}} nodes, {.val {x$n_periods}} time periods ({net_type})")
	}

	# model config
	config_parts <- c(x$family, x$method)
	if (isTRUE(x$fix_receiver) && !isTRUE(x$symmetric)) {
	  config_parts <- c(config_parts, "fix_receiver")
	}
	cli::cli_text("Config: {.field {paste(config_parts, collapse = ' | ')}}")

	# status line
	n_obs <- if (!is.null(x$nobs)) x$nobs else NA
	if (x$convergence) {
	cli::cli_text("Status: {.field converged} | N = {.val {n_obs}} | Log-Lik: {.val {round(x$ll, 2)}} | AIC: {.val {round(AIC(x), 1)}}")
	} else {
	cli::cli_text("Status: {.emph not converged} | N = {.val {n_obs}} | Log-Lik: {.val {round(x$ll, 2)}}")
	}

	# coefficients as formatted table
	if (nrow(x$summ) > 0) {
	cli::cli_text("\nCoefficients:")
	have_se <- "se" %in% colnames(x$summ) && any(is.finite(x$summ$se))
	if (have_se) {
	  coef_mat <- as.matrix(x$summ[, c("coef", "se")])
	  colnames(coef_mat) <- c("Estimate", "Std. Err")
	  print(round(coef_mat, digits))
	} else {
	  coef_mat <- as.matrix(x$summ[, "coef", drop = FALSE])
	  colnames(coef_mat) <- "Estimate"
	  print(round(coef_mat, digits))
	  if (isTRUE(x$full_bilinear)) {
		  cli::cli_text("{.emph (SEs not computed; use boot_sir(type = \"dyad\") for full-bilinear bipartite fits)}")
	  } else if (isFALSE(x$se_reliable)) {
		  cli::cli_text("{.emph (SEs unavailable or unreliable; use summary() for diagnostics or boot_sir() for resampling)}")
	  } else {
		  cli::cli_text("{.emph (SEs not computed)}")
	  }
	}
	}

	cli::cli_text("\nUse {.code summary()} for detailed results")
	invisible(x)
}

# check that newdata components share the time dimension implied by X and
# abort on a mismatch; internal, not exported
validate_predict_newdata <- function(object, W, X, Z, W_recv = NULL) {
	n1 <- if (!is.null(object$n1)) object$n1 else dim(object$Y)[1]
	n2 <- if (!is.null(object$n2)) object$n2 else dim(object$Y)[2]
	p_fit <- if (!is.null(object$p)) object$p else 0L
	p2_fit <- if (!is.null(object$p2)) object$p2 else 0L
	q_fit <- if (!is.null(object$q)) object$q else length(object$theta)

	if (is.null(X)) {
		cli::cli_abort(c(
			"{.field newdata} must supply (or the fit must store) the predictor array {.field X}.",
			"i" = "Provide {.code newdata = list(W = , X = , Z = )}."
		))
	}
	if (!is.array(X) || length(dim(X)) != 3L) {
		cli::cli_abort("{.field newdata$X} must be a 3D array.")
	}
	if (dim(X)[1] != n1 || dim(X)[2] != n2) {
		cli::cli_abort(c(
			"{.field newdata$X} has incompatible first two dimensions.",
			"x" = "Expected {.val {n1}} x {.val {n2}}, got {.val {dim(X)[1]}} x {.val {dim(X)[2]}}."
		))
	}

		if (p_fit > 0L) {
			if (is.null(W)) {
				cli::cli_abort("{.field newdata$W} is required because the fit estimated influence covariates.")
			}
			if (!is.array(W) || !(length(dim(W)) %in% c(3L, 4L))) {
				cli::cli_abort("{.field newdata$W} must be a 3D or 4D array.")
			}
			if (!is.null(object$W_recv) && length(dim(W)) != 3L) {
				cli::cli_abort(c(
					"{.field newdata$W} must be a 3D static sender-side array for full-bilinear bipartite fits.",
					"i" = "Dynamic full-bilinear bipartite prediction is not implemented; fit separate static scenarios or omit {.field newdata$W} to reuse the fitted array."
				))
			}
		if (dim(W)[1] != n1 || dim(W)[2] != n1) {
			cli::cli_abort(c(
				"{.field newdata$W} has incompatible first two dimensions.",
				"x" = "Expected sender/influence dimensions {.val {n1}} x {.val {n1}}, got {.val {dim(W)[1]}} x {.val {dim(W)[2]}}."
			))
		}
		if (dim(W)[3] != p_fit) {
			cli::cli_abort(c(
				"{.field newdata$W} carries {dim(W)[3]} influence covariate{?s} but the fit estimated {p_fit}.",
				"i" = "Supply W with the same covariate dimension as the fitted model."
			))
		}
		} else if (!is.null(W)) {
			cli::cli_abort("{.field newdata$W} was supplied, but the fit has no influence-covariate terms.")
		}

	if (!is.null(object$W_recv)) {
		if (is.null(W_recv)) {
			cli::cli_abort("{.field newdata$W_recv} is required for a full-bilinear bipartite fit.")
		}
		if (!is.array(W_recv) || length(dim(W_recv)) != 3L) {
			cli::cli_abort("{.field newdata$W_recv} must be a 3D array.")
		}
		if (dim(W_recv)[1] != n2 || dim(W_recv)[2] != n2) {
			cli::cli_abort(c(
				"{.field newdata$W_recv} has incompatible first two dimensions.",
				"x" = "Expected receiver/influence dimensions {.val {n2}} x {.val {n2}}, got {.val {dim(W_recv)[1]}} x {.val {dim(W_recv)[2]}}."
			))
		}
		if (dim(W_recv)[3] != p2_fit) {
			cli::cli_abort(c(
				"{.field newdata$W_recv} carries {dim(W_recv)[3]} receiver influence covariate{?s} but the fit estimated {p2_fit}.",
				"i" = "Supply W_recv with the same covariate dimension as the fitted model."
			))
		}
	} else if (!is.null(W_recv)) {
		cli::cli_abort("{.field newdata$W_recv} was supplied, but the fit has no receiver-side influence covariates.")
	}

	if (q_fit > 0L) {
		if (is.null(Z)) {
			cli::cli_abort("{.field newdata$Z} is required because the fit estimated direct covariates.")
		}
		if (!is.array(Z) || !(length(dim(Z)) %in% c(3L, 4L))) {
			cli::cli_abort("{.field newdata$Z} must be a 3D or 4D array.")
		}
		if (dim(Z)[1] != n1 || dim(Z)[2] != n2) {
			cli::cli_abort(c(
				"{.field newdata$Z} has incompatible first two dimensions.",
				"x" = "Expected {.val {n1}} x {.val {n2}}, got {.val {dim(Z)[1]}} x {.val {dim(Z)[2]}}."
			))
		}
	} else if (!is.null(Z)) {
		cli::cli_abort("{.field newdata$Z} was supplied, but the fit has no direct covariate terms.")
	}

	# the bilinear computation takes its time span from X (m x m x T)
	n_periods <- if (length(dim(X)) >= 3) dim(X)[3] else 1L
	t_dims <- c(X = n_periods)
	# Z is m x m x q x T (4D) or m x m x T (3D, q = 1); trailing dim is time
	if (!is.null(Z)) {
		z_t <- dim(Z)[length(dim(Z))]
		t_dims <- c(t_dims, Z = z_t)
	}
	# dynamic (time-varying) W is a 4D array; trailing dim is time
	if (!is.null(W) && length(dim(W)) == 4) {
		t_dims <- c(t_dims, W = dim(W)[4])
	}
	if (length(unique(t_dims)) > 1L) {
		detail <- paste(names(t_dims), t_dims, sep = " = ", collapse = ", ")
		cli::cli_abort(c(
			"Time dimensions in {.field newdata} disagree.",
			"x" = "Component time lengths: {detail}.",
			"i" = "All time-varying inputs must span the same number of periods; values are not recycled."
		))
	}
	# theta (covariate) dimension lives in the 3rd slot of a 4D Z (m x m x q x T)
	if (!is.null(Z)) {
		q_new <- if (length(dim(Z)) >= 4) dim(Z)[3] else 1L
		if (q_new != q_fit) {
			cli::cli_abort(c(
				"{.field newdata$Z} carries {q_new} covariate{?s} but the fit estimated {q_fit} {.field theta} coefficient{?s}.",
				"i" = "Supply a {.field Z} whose covariate dimension matches the fitted covariates."
			))
		}
	}
	invisible(TRUE)
}

.sir_prediction_mask <- function(arr, object) {
	if (isTRUE(object$bipartite) || object$n1 != object$n2) return(arr)
	if (isTRUE(object$symmetric)) {
		for (tt in seq_len(dim(arr)[3])) {
			upper <- arr[, , tt]
			upper[lower.tri(upper, diag = TRUE)] <- NA
			upper0 <- upper
			upper0[is.na(upper0)] <- 0
			sym <- upper0 + t(upper0)
			diag(sym) <- NA
			arr[, , tt] <- sym
		}
		return(arr)
	}
	set_square_diagonal(arr, NA)
}

#' Predictions from a Fitted SIR Model
#'
#' Generates predictions from a fitted SIR model for the training data or
#' for new data. Predictions can be on the link scale (linear predictor) or
#' the response scale (expected counts, probabilities, or means).
#'
#' For model-implied scenario analysis, supply modified arrays in
#' \code{newdata}. For example, to see how fitted values change when a covariate
#' is increased by one unit, pass the modified Z array while keeping W and X
#' from the original fit. Causal counterfactual interpretation requires
#' additional design assumptions.
#'
#' @param object A fitted \code{sir} object from \code{\link{sir}}.
#' @param newdata Optional named list with components \code{W} (3D or 4D
#'   array; 3D only for full-bilinear bipartite fits), \code{X} (3D array),
#'   \code{Z} (3D or 4D array), and for
#'   full-bilinear bipartite fits \code{W_recv} (3D receiver-side array), for
#'   scenario prediction. Dimensions must match the original fit. Any
#'   component not supplied is taken from the original fit. All supplied
#'   time-varying components must agree on the number of time periods, and
#'   \code{Z} must carry the same number of covariates as the fit; a genuine
#'   mismatch is an error rather than being silently recycled. If NULL
#'   (default), returns predictions for the training data. Note: unlike many R
#'   predict methods, \code{newdata} is a list of arrays, not a data frame.
#'   For a full-bilinear bipartite fit, supply \code{newdata$W_recv} to vary the
#'   receiver-side influence structure; otherwise the fitted \code{W_recv} is
#'   reused. Unlike the square one-mode case the bipartite diagonal is a genuine
#'   prediction and is not set to NA.
#' @param type Character string: \code{"link"} for linear predictor (eta)
#'   or \code{"response"} for expected values on the original scale.
#'   Default is \code{"response"}.
#' @param ... Additional arguments (unused).
#' @return An array (n1 x n2 x T) of predicted values on the requested scale.
#'
#' @examples
#' dat <- sim_sir(m = 10, T_len = 20, p = 2, q = 1, family = "poisson", seed = 1)
#' fit <- sir(dat$Y, W = dat$W, X = dat$X, Z = dat$Z, family = "poisson", seed = 1)
#' # in-sample fitted values (response scale)
#' pred <- predict(fit)
#' # scenario: only W/X/Z are read from newdata (Y is ignored)
#' Zcf <- dat$Z; Zcf[, , 1, ] <- Zcf[, , 1, ] + 1
#' pred_cf <- predict(fit, newdata = list(W = dat$W, X = dat$X, Z = Zcf))
#' @export
predict.sir <- function(object, newdata = NULL,
						type = c("response", "link"), ...) {
	type <- match.arg(type)

	# full-bilinear bipartite fits carry a receiver-side W_recv and need the
	# two-sided linear predictor; B is not the identity here.
	if (!is.null(object$W_recv)) {
		if (!is.null(newdata)) {
			W_new <- if (!is.null(newdata$W)) newdata$W else object$W
			X_new <- if (!is.null(newdata$X)) newdata$X else object$X
			Z_new <- if (!is.null(newdata$Z)) newdata$Z else object$Z
			W_recv_new <- if (!is.null(newdata$W_recv)) newdata$W_recv else object$W_recv
			validate_predict_newdata(object, W_new, X_new, Z_new, W_recv_new)
			eta <- eta_tab_bipartite(object$tab, W_new, W_recv_new,
									 X_new, Z_new, object$p, object$p2, object$q)
		} else {
			if (type == "response") return(object$fitted.values)
			eta <- eta_tab_bipartite(object$tab, object$W, object$W_recv,
									 object$X, object$Z, object$p, object$p2, object$q)
		}
		if (type == "response") {
			return(switch(object$family,
				poisson = exp(eta),
				binomial = 1 / (1 + exp(-eta)),
				eta))
		}
		return(eta)
	}

	fr <- isTRUE(object$fix_receiver)

	if (!is.null(newdata)) {
	# extract components from newdata, falling back to original fit. any
		# component omitted from newdata is taken from the fit
			W_new <- if (!is.null(newdata$W)) newdata$W else if (object$p > 0L) object$W else NULL
		if (!is.null(W_new) && !isTRUE(object$bipartite) &&
			object$n1 == object$n2) {
			W_new <- set_square_diagonal(W_new, 0)
		}
		X_new <- if (!is.null(newdata$X)) newdata$X else object$X
		if (!is.null(X_new) && !isTRUE(object$bipartite) &&
			object$n1 == object$n2) {
			X_new <- set_square_diagonal(X_new, 0)
		}
		Z_new <- if (!is.null(newdata$Z)) newdata$Z else object$Z

		# guard against silent time-dimension recycling before computing eta
			validate_predict_newdata(object, W_new, X_new, Z_new)

	# calculate linear predictor
	eta <- eta_tab(object$tab, W_new, X_new, Z_new, fix_receiver=fr)

		if (type == "response") {
		  if (object$family == "poisson") {
			out <- exp(eta)
		  } else if (object$family == "binomial") {
			out <- 1 / (1 + exp(-eta))
		  } else {
			out <- eta
		  }
			  return(.sir_prediction_mask(out, object))
			} else {
			  return(.sir_prediction_mask(eta, object))
			}
		} else {
		# return predictions from training data
		if (type == "response") {
			return(object$fitted.values)
		} else {
				eta <- eta_tab(object$tab, object$W, object$X, object$Z, fix_receiver=fr)
				return(.sir_prediction_mask(eta, object))
			}
		}
}

#' Variance-Covariance Matrix for SIR Model Parameters
#'
#' Returns the variance-covariance matrix of the estimated parameters. Several
#' types are available; for relational (dyadic) data the cluster-robust types are
#' usually more cautious because the classical and HC0-robust covariances assume
#' independent dyad-periods and often undercover.
#'
#' @param object A fitted \code{sir} object from \code{\link{sir}}.
#' @param type Character string:
#'   \itemize{
#'     \item \code{"cluster"} (default) --- multiway cluster-robust covariance
#'       on sender, receiver, and time margins for directed-network data when
#'       the Hessian bread is stable.
#'     \item \code{"classical"} --- inverse-Hessian covariance.
#'     \item \code{"robust"} --- HC0 sandwich; corrects heteroskedasticity /
#'       overdispersion only, \emph{not} dyadic dependence.
#'     \item \code{"twoway"} --- alias for \code{"cluster"}.
#'     \item \code{"dyad"} --- clusters the directed dyad across time.
#'   }
#'   The cluster types require the classical covariance as their bread
#'   (\code{calc_se = TRUE}, the default) and are unavailable for dynamic (4D)
#'   \code{W} and for full-bilinear bipartite fits (use
#'   \code{boot_sir(type = "dyad")} there).
#' @param ... Additional arguments (unused).
#' @return A square matrix with rows and columns named by parameter.
#'   Returns NULL if standard errors were not computed (\code{calc_se = FALSE}).
#' @seealso \code{\link{confint.sir}} for confidence intervals,
#'   \code{\link{boot_sir}} for resampling-based inference.
#' @export
	vcov.sir <- function(object, type = c("cluster", "classical", "robust", "twoway", "dyad"), ...) {
		type <- match.arg(type)
		if (!isTRUE(object$convergence)) {
			cli::cli_abort(c(
				"Model did not converge; Wald covariance estimates are not reliable.",
				"i" = "Refit until convergence before reporting standard errors or intervals."
			))
		}
		# full-bilinear bipartite fits have no analytic covariance; the only valid
	# inference is the dyad jackknife.
	if (isTRUE(object$full_bilinear)) {
		cli::cli_abort(c(
			"Analytic covariance is not available for full-bilinear bipartite fits.",
			"i" = "Use {.code boot_sir(fit, type = \"dyad\")} and read SEs/intervals from that result."
		))
	}
	if (isFALSE(object$se_reliable) && type %in% c("classical", "robust")) {
		cli::cli_abort(c(
			"Requested {.val {type}} covariance is based on an unreliable Hessian for this fit.",
			"i" = "Use {.code boot_sir()} or refit/simplify the model before reporting Wald inference."
		))
	}
			if (type %in% c("cluster", "twoway", "dyad")) {
		  # cluster-robust sandwich for dyadic dependence
	  by <- if (type == "dyad") "dyad" else "twoway"
	  return(.sir_vcov_cluster(object, by = by))
	}
	if (type == "robust") {
	  V <- object$vcov_robust
	  if (is.null(V)) {
		  cli::cli_abort(c(
			  "Robust covariance matrix not available.",
			  "i" = "Refit with {.code calc_se = TRUE}; robust SEs are unavailable for some paths (e.g. when the sandwich could not be formed)."
		  ))
	  }
	  return(V)
	}
	if (is.null(object$vcov)) {
	  cli::cli_abort(c(
		  "Classical covariance matrix not available.",
		  "i" = "Refit with {.code calc_se = TRUE}."
	  ))
	}
	object$vcov
}

#' Confidence Intervals for SIR Model Parameters
#'
#' Computes confidence intervals using either Wald-based intervals (cluster-
#' robust by default for supported static fits) or bootstrap percentile
#' intervals. Wald intervals require that the model was fit with
#' \code{calc_se = TRUE}.
#' Bootstrap intervals require a \code{\link{boot_sir}} result and tend to
#' be more reliable when the Hessian is ill-conditioned.
#'
#' @param object A fitted \code{sir} object from \code{\link{sir}}.
#' @param parm Character vector of parameter names or integer indices to
#'   include. If NULL (default), returns intervals for all parameters.
#' @param level Confidence level between 0 and 1. Default is 0.95.
#' @param boot Optional \code{boot_sir} object from \code{\link{boot_sir}}.
#'   When provided, percentile intervals from the bootstrap distribution are
#'   used instead of Wald intervals.
#' @param se.type Character. For Wald intervals, which standard errors to use:
#'   \code{"cluster"} (default; multiway sender, receiver, and time clustering
#'   for supported static fits), \code{"classical"} (inverse-Hessian),
#'   \code{"robust"} (HC0 sandwich), \code{"twoway"} (alias for
#'   \code{"cluster"}), or \code{"dyad"} (directed-dyad clustering). Cluster
#'   intervals still rely on the Hessian bread being stable. Ignored
#'   when \code{boot} is supplied. Non-classical types require \code{calc_se =
#'   TRUE}; cluster types are unavailable for dynamic (4D) \code{W}.
#' @param ... Additional arguments (unused).
#' @return A matrix with one row per parameter and columns for the lower
#'   and upper bounds, labeled by percentage (e.g., \code{"2.5 \%"} and
#'   \code{"97.5 \%"}).
#' @seealso \code{\link{boot_sir}} for bootstrap inference,
#'   \code{\link{vcov.sir}} for the variance-covariance matrix.
#' @export
	confint.sir <- function(object, parm = NULL, level = 0.95, boot = NULL,
							se.type = c("cluster", "classical", "robust", "twoway", "dyad"), ...) {
		se.type <- match.arg(se.type)
		if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
			level <= 0 || level >= 1) {
			cli::cli_abort("{.arg level} must be a single number between 0 and 1.")
		}
		cf <- coef(object)
	pnames <- names(cf)
	a <- (1 - level) / 2
	pct <- paste0(format(100 * c(a, 1 - a), trim = TRUE, digits = 3), " %")

	# full-bilinear bipartite fits have no analytic SEs; require a boot result.
	if (is.null(boot) && isTRUE(object$full_bilinear)) {
		cli::cli_abort(c(
			"Wald intervals are not available for full-bilinear bipartite fits (no analytic SEs).",
			"i" = "Run {.code b <- boot_sir(fit, type = \"dyad\")} and pass {.code confint(fit, boot = b)}."
		))
	}

	if (!is.null(boot)) {
		if (!inherits(boot, "boot_sir")) {
			cli::cli_abort("{.arg boot} must be a {.cls boot_sir} object from {.fn boot_sir}.")
		}
		if (!identical(boot$family, object$family)) {
			cli::cli_abort("{.arg boot} was computed for family {.val {boot$family}}, but the model uses {.val {object$family}}.")
		}
		if (!identical(boot$param_names, pnames)) {
			cli::cli_abort("{.arg boot} parameter names do not match this model.")
		}
		if (!isTRUE(all.equal(unname(boot$point_est), unname(cf), tolerance = 1e-8))) {
			cli::cli_abort("{.arg boot} was computed from a different fitted model (point estimates differ).")
		}
		# bootstrap percentile / jackknife intervals
		ci <- confint(boot, level = level)[, , drop = FALSE]
	ci <- ci[pnames, , drop = FALSE]
	rownames(ci) <- pnames
	colnames(ci) <- pct
	} else {
		# wald intervals; pick classical, robust (HC0), or cluster-robust SEs.
		# use cluster-robust standard errors
		if (!isTRUE(object$convergence)) {
		  cli::cli_abort(c(
			"Model did not converge; Wald intervals are not reliable.",
			"i" = "Refit until convergence, or use a clearly labeled resampling sensitivity check."
		  ))
		}
		if (isFALSE(object$se_reliable) && se.type %in% c("classical", "robust")) {
	  cli::cli_abort(c(
		"Requested {.val {se.type}} Wald intervals are based on an unreliable Hessian for this fit.",
		"i" = "Use {.code boot_sir()} or refit/simplify the model before reporting Wald inference."
	  ))
	}
	if (se.type == "classical") {
	  se <- object$summ$se
	  if (is.null(se) || all(is.na(se))) {
		cli::cli_abort("No classical standard errors available. Refit with {.code calc_se = TRUE} or provide bootstrap results via the {.arg boot} argument.")
	  }
	} else if (se.type == "robust") {
	  se <- object$summ$rse
	  if (is.null(se) || all(is.na(se))) {
		cli::cli_abort(c(
		  "No robust standard errors available for robust Wald intervals.",
		  "i" = "Refit with {.code calc_se = TRUE}, or use {.code se.type = \"classical\"} or the {.arg boot} argument."
		))
	  }
	} else {
	  by <- if (se.type == "dyad") "dyad" else "twoway"
	  se <- sqrt(diag(.sir_vcov_cluster(object, by = by)))
	}
	z <- qnorm(1 - a)
	ci <- cbind(cf - z * se, cf + z * se)
	rownames(ci) <- pnames
	colnames(ci) <- pct
	}

	if (!is.null(parm)) {
	if (is.character(parm)) {
	  bad <- setdiff(parm, pnames)
	  if (length(bad)) cli::cli_abort("Unknown parameter{?s} in {.arg parm}: {.val {bad}}.")
	} else {
	  if (any(parm < 1 | parm > length(cf))) {
		cli::cli_abort("{.arg parm} index out of range (model has {length(cf)} parameter{?s}).")
	  }
	}
	ci <- ci[parm, , drop = FALSE]
	}

	ci
}
